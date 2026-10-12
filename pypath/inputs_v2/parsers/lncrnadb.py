"""
Read the lncRNAdb association export.

Source file (XML), the response of the lncRNAdb REST service for all entries:

    lncrnadb_association.xml  (``<Search_Results>`` of ``<Results>``)

Each ``<Results>`` element is one lncRNA: ``<nomenclature>`` (``Name``, ``Alias``,
``Biotype``), ``<species>`` (one ``<Entry Species="...">`` per species in which the
lncRNA is recorded, the text being a genome locus or ``N/A``) and ``<association>``,
which holds the ``<association>`` elements: the interaction partner with its
``ComponentType`` and ``ComponentID``, a ``description`` and a ``pubmedID``. An
association is not tied to one of the species of its lncRNA. An entry without
associations holds an empty ``<association/>``; an association without a partner
(empty type and ID, PubMed ID ``None``) is a placeholder.
"""

from __future__ import annotations

import re
from collections.abc import Generator
from typing import Any

from pypath.inputs_v2._publications import is_pmid
from pypath.inputs_v2.parsers.base import read_xml_root

# The source writes a missing value in several ways.
_NULLS = frozenset({'', 'None', 'NULL', 'N/A', 'None known'})
# A PubMed ID is a number of seven or eight digits. The PubMed ID cell
# sometimes holds a status word instead (``Active``, ``Validated``,
# ``Confirmed``), two IDs, or an ID after a colon.


def _text(element: Any, path: str) -> str | None:
    value = (element.findtext(path) or '').strip()
    return None if value in _NULLS else value


def iter_associations(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield the associations of every lncRNA entry, with the entry's fields.

    Fields: ``lncrna`` (``Name``), ``aliases`` (the ``Alias`` cell split at commas
    and semicolons), ``biotype``, ``species`` (the ``Species`` attributes, as
    written, one per ``<Entry>``), ``partner_type`` (``ComponentType``, as
    written), ``partner`` (``ComponentID``), ``description`` and ``pmids`` (the
    PubMed IDs of the ``pubmedID`` cell, which holds one or two numbers or a status word). Entries without
    associations yield nothing.
    """
    root = read_xml_root(opener)
    if root is None or root.tag != 'Search_Results':
        raise ValueError('lncRNAdb download is not a Search_Results document')
    entries = 0
    count = 0
    for entry in root.findall('Results'):
        entries += 1
        nomenclature = entry.find('nomenclature')
        name = _text(nomenclature, 'Name') if nomenclature is not None else None
        if name is None:
            raise ValueError('lncRNAdb entry without a name')
        aliases = _text(nomenclature, 'Alias')
        species = [
            (e.get('Species') or '').strip() for e in entry.findall('species/Entry')
        ]
        for association in entry.findall('association/association'):
            count += 1
            yield {
                'lncrna': name,
                'aliases': [a.strip() for a in re.split(r'[,;]', aliases) if a.strip()] if aliases else [],
                'biotype': _text(nomenclature, 'Biotype'),
                'species': [s for s in species if s],
                'partner_type': _text(association, 'ComponentType'),
                'partner': _text(association, 'ComponentID'),
                'description': _text(association, 'description'),
                'pmids': _pmid_list(association.findtext('pubmedID')),
            }
    if not entries:
        raise ValueError('lncRNAdb download has no entries')
    if not count:
        raise ValueError('lncRNAdb download has no associations')


def _pmid_list(value: str | None) -> list[str]:
    return [token for token in re.split(r'[\s,;:]+', value or '') if is_pmid(token)]
