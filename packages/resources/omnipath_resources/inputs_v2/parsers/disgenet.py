"""
Readers for the DisGeNET gene-disease association file.

``curated_gene_disease_associations.tsv`` is a gzipped tab separated table,
one gene-disease association per line, with a header row. Numbers are padded
with blanks (``geneId``, ``NofPmids``). The shared TSV reader reads it; the
functions here check the header, trim the cells and give the lists.
"""

from __future__ import annotations

from collections.abc import Generator
from typing import Any

from pypath.inputs_v2.parsers.base import iter_tsv

__all__ = ['iter_associations', 'iter_diseases']

COLUMNS = (
    'geneId',
    'geneSymbol',
    'DSI',
    'DPI',
    'diseaseId',
    'diseaseName',
    'diseaseType',
    'diseaseClass',
    'diseaseSemanticType',
    'score',
    'EI',
    'YearInitial',
    'YearFinal',
    'NofPmids',
    'NofSnps',
    'source',
)


def _cell(value: str | None) -> str | None:
    value = None if value is None else value.strip()
    return value or None


def iter_associations(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield one record per gene-disease association of the file.

    Output fields, the columns of the file with blanks trimmed and empty cells
    as None:

    - ``entrez_id``, ``symbol``: the gene (``geneId``, ``geneSymbol``).
    - ``disease_id``, ``disease_name``, ``disease_type``: the disease (a UMLS
      concept unique identifier such as ``C0036341``), its name, and its type
      as DisGeNET gives it (``disease``, ``phenotype``, ``group``).
    - ``dsi``, ``dpi``: the disease specificity and pleiotropy indices of the
      gene.
    - ``score``: the DisGeNET score of the association (0 to 1); ``ei``: its
      evidence index.
    - ``year_initial``, ``year_final``: the years of the first and the last
      publication that report it.
    - ``n_pmids``, ``n_snps``: the number of publications and of variants that
      support it.
    - ``sources``: the list of the databases that report it (``CTD_human``,
      ``UNIPROT``).
    """
    for row in iter_tsv(opener, required=COLUMNS):
        sources = _cell(row['source'])
        yield {
            'entrez_id': _cell(row['geneId']),
            'symbol': _cell(row['geneSymbol']),
            'disease_id': _cell(row['diseaseId']),
            'disease_name': _cell(row['diseaseName']),
            'disease_type': _cell(row['diseaseType']),
            'dsi': _cell(row['DSI']),
            'dpi': _cell(row['DPI']),
            'score': _cell(row['score']),
            'ei': _cell(row['EI']),
            'year_initial': _cell(row['YearInitial']),
            'year_final': _cell(row['YearFinal']),
            'n_pmids': _cell(row['NofPmids']),
            'n_snps': _cell(row['NofSnps']),
            'sources': sources.split(';') if sources else [],
        }


def iter_diseases(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield each disease of the file once, in order of first appearance.

    Output fields: ``disease_id``, ``name`` and ``type`` as in
    :func:`iter_associations`. A disease with two names or two types raises.
    """
    seen: dict[str, tuple[str | None, str | None]] = {}
    for row in iter_tsv(opener, required=COLUMNS):
        disease = _cell(row['diseaseId'])
        value = (_cell(row['diseaseName']), _cell(row['diseaseType']))
        if disease is None or value[0] is None:
            raise ValueError(f'DisGeNET association without a disease or a name: {row!r}')
        if disease not in seen:
            seen[disease] = value
            yield {'disease_id': disease, 'name': value[0], 'type': value[1]}
        elif seen[disease] != value:
            raise ValueError(f'DisGeNET disease {disease} has two names or types: {seen[disease]!r}, {value!r}')
