"""
Read the interaction table of LncRNADisease.

Source file (tab separated, no header row, no quoting):

    data2.txt (LncRNADisease_data2.txt)

A row is one literature-curated interaction of a lncRNA with another molecule, with
ten columns: the row number, the lncRNA name, the name of the interacting molecule,
the interaction level, the interaction class, a sentence from the paper, the species,
the aliases of the lncRNA (``; `` separated), a sequence accession of the lncRNA and
the PubMed ID. The level names the molecule types of the two partners, as
``<first>-<second>``: ``RNA-Protein``, ``RNA-RNA``, ``RNA-DNA``, ``DNA-TF``,
``DNA-DNA``, ``DNA-Protein``. A first partner of ``DNA`` is the gene of the lncRNA (the
paper describes the ``DNA-TF`` level as the promoters of lncRNA genes binding and
being regulated by transcription factors). ``N/A`` marks a missing value.
"""

from __future__ import annotations

from collections.abc import Generator
from typing import Any

from pypath.inputs_v2._publications import is_pmid
from pypath.inputs_v2.parsers.base import iter_tsv

# Columns of the table, in source order.
COLUMNS = (
    'id',
    'lncrna',
    'partner',
    'level',
    'interaction_class',
    'description',
    'species',
    'aliases',
    'accession',
    'pmid',
)

_NULLS = frozenset({'', 'N/A'})


def _cell(value: Any) -> str | None:
    value = (value or '').strip()
    return None if value in _NULLS else value


def iter_interactions(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield the rows of the interaction table with their cells split and cleaned.

    Fields: ``id``, ``lncrna`` and ``partner`` (names), ``level`` (as in the
    source) with ``first_type`` and ``second_type`` (its two parts; None for
    ``N/A``), ``interaction_class``, ``description``, ``species`` (as in the
    source), ``aliases`` (a list), ``accession`` and ``pmids`` (a list of PubMed
    IDs).
    """
    count = 0
    for row in iter_tsv(opener, quotechar=None, fieldnames=COLUMNS, strict=True):
        count += 1
        level = _cell(row['level'])
        first, _, second = (level or '').partition('-')
        if level is not None and not (first and second):
            raise ValueError(f'LncRNADisease interaction level is not <first>-<second>: {level!r}')
        pmid = _cell(row['pmid'])
        if pmid is not None and not pmid.isdigit():
            raise ValueError(f'LncRNADisease PMID is not a number: {pmid!r}')
        aliases = _cell(row['aliases'])
        yield {
            'id': _cell(row['id']),
            'lncrna': _cell(row['lncrna']),
            'partner': _cell(row['partner']),
            'level': level,
            'first_type': first or None,
            'second_type': second or None,
            'interaction_class': _cell(row['interaction_class']),
            'description': _cell(row['description']),
            'species': _cell(row['species']),
            'aliases': [a.strip() for a in aliases.split(';') if a.strip()] if aliases else [],
            'accession': _cell(row['accession']),
            'pmids': [pmid] if is_pmid(pmid) else [],
        }
    if not count:
        raise ValueError('LncRNADisease interaction table has no rows')
