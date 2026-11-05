"""Transporter-substrate pairs of the substrate-based ontology of human solute carriers.

Meixner et al. 2020 (Mol Syst Biol 16: e9652) annotate 446 human SLCs (Table EV1,
sheet ``slc_annotation``) with their substrates and coupled ions as names, and map
each name to ChEBI (Table EV2, sheet ``chebi_mapping``). A pair is a transporter
row and one name of its ``Substrates`` or ``Coupled_ions`` cell, as the legacy
``slc_interactions`` read them; names the mapping gives no ChEBI identifier for are
left out, as legacy left them out.
"""

from __future__ import annotations

from collections.abc import Generator
import re
from typing import Any

from pypath.inputs_v2.parsers.base import iter_excel

ANNOTATION_SHEET = 'slc_annotation'
MAPPING_SHEET = 'chebi_mapping'
# Role, the column of names, the column of their references.
ROLES = (
    ('substrate', 'Substrates', 'Substrates_PMID'),
    ('coupled_ion', 'Coupled_ions', 'Coupled_ions_PMID'),
)
# Cells the tables use for "no value".
NOTHING = frozenset({'', 'Unknown', 'None'})
# Reference cells mix PMIDs with notes ("12133831 (rabbit)", "11827462; 18337592").
PMID = re.compile(r'(?<!\d)\d{5,9}(?!\d)')


def _value(value: Any) -> str | None:
    text = ' '.join(str(value).split()) if value is not None else ''
    return None if text in NOTHING else text


def pmids(value: Any) -> list[str]:
    """The PubMed identifiers of a reference cell, without its notes."""
    return sorted(set(PMID.findall(str(value or ''))))


def chebi_mapping(opener) -> dict[str, dict[str, str]]:
    """Substrate name -> ChEBI identifier and term, for names that have one."""
    mapping = {}
    for row in iter_excel(
        opener,
        sheet=MAPPING_SHEET,
        required=('Substrate_name_annotation', 'chebi_id', 'chebi_term'),
    ):
        name, chebi = row.get('Substrate_name_annotation'), _value(row.get('chebi_id'))
        if name and chebi:
            mapping[name] = {'chebi': chebi, 'chebi_term': _value(row.get('chebi_term'))}
    return mapping


def transport_substrates(
    opener,
    *,
    mapping_download,
    **kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """One record per transporter and substrate or coupled ion.

    ``mapping_download`` is the Table EV2 download, injected by the ``slc`` module;
    ``force_refresh`` is forwarded to it.
    """
    if not opener or not opener.result:
        return
    mapping = chebi_mapping(mapping_download.open(force_refresh=kwargs.get('force_refresh', False)))
    for row in iter_excel(
        opener,
        sheet=ANNOTATION_SHEET,
        required=('HGNC_symbol', *(column for _, names, refs in ROLES for column in (names, refs))),
    ):
        if not (genesymbol := _value(row.get('HGNC_symbol'))):
            continue
        for role, names_column, refs_column in ROLES:
            cell = row.get(names_column)
            if cell is None:
                continue
            for name in str(cell).split('; '):
                if name in NOTHING or not (substrate := mapping.get(name)):
                    continue
                yield {
                    'genesymbol': genesymbol,
                    'entrez': _value(row.get('Entrez_ID')),
                    'ensg': _value(row.get('Ensembl_ID')),
                    'role': role,
                    'substrate_name': name,
                    'chebi': substrate['chebi'],
                    'pmids': pmids(row.get(refs_column)),
                    'substrate_class': _value(row.get('Substrate_class')),
                    'transport_mechanism': _value(row.get('Transport_mechanism')),
                    'localization': _value(row.get('Subcellular_localization')),
                }
