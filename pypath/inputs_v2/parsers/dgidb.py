"""
Readers for the DGIdb downloads (``categories.tsv``, ``genes.tsv``, ``interactions.tsv``).

All three files are tab separated, with two comment lines first (``# Data
version: 2026-09`` and ``# DGIdb version: v.5.0.13``) and then a header row.
The shared TSV reader reads them; the functions here are the DGIdb
interpretation: they check the header, give empty cells as ``None``, join the
gene categories to the gene concepts of ``genes.tsv`` and give one record per
interaction type of a drug-gene interaction claim.
"""

from __future__ import annotations

from collections.abc import Generator
import functools
from typing import Any

from pypath.inputs_v2.base import Download
from pypath.inputs_v2.parsers.base import iter_tsv

__all__ = ['iter_categories', 'iter_category_terms', 'iter_interactions']

# The two comment lines before the header row.
PREAMBLE_LINES = 2
CATEGORY_COLUMNS = (
    'gene_claim_name',
    'gene_category_name',
    'source_db_name',
    'source_db_version',
)
GENE_COLUMNS = (
    'gene_claim_name',
    'nomenclature',
    'concept_id',
    'gene_name',
    'source_db_name',
    'source_db_version',
)
INTERACTION_COLUMNS = (
    'gene_claim_name',
    'gene_concept_id',
    'gene_name',
    'drug_claim_name',
    'drug_concept_id',
    'drug_name',
    'drug_is_approved',
    'drug_is_immunotherapy',
    'drug_is_antineoplastic',
    'interaction_source_db_name',
    'interaction_source_db_version',
    'interaction_types',
    'interaction_score',
    'drug_specificity_score',
    'gene_specificity_score',
    'evidence_score',
)

# The files hold quotes as data (``""`` stands for an empty list of types), so
# the reader does not interpret them.
_read = functools.partial(iter_tsv, skip_rows=PREAMBLE_LINES, quotechar=None)


def _cell(value: str | None) -> str | None:
    """A cell without surrounding blanks; an empty cell, ``""`` or a missing one is None."""
    if value is None:
        return None
    value = value.strip()
    return None if value in ('', '""') else value


def iter_category_terms(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield each gene category of ``categories.tsv`` once, in order of first appearance.

    Output fields: ``category``, the category as the file writes it
    (``KINASE``, ``DRUGGABLE GENOME``).
    """
    seen: set[str] = set()
    for row in _read(opener, required=CATEGORY_COLUMNS):
        category = _cell(row['gene_category_name'])
        if category is None:
            raise ValueError(f'DGIdb categories: a claim without category: {row!r}')
        if category not in seen:
            seen.add(category)
            yield {'category': category}


def _gene_concepts(genes_download: Download, **kwargs: Any) -> dict[tuple[str, str], tuple[str, str]]:
    """(gene claim name, source) -> (concept id, gene name), from ``genes.tsv``.

    A claim of a source that maps to two different concepts raises.
    """
    concepts: dict[tuple[str, str], tuple[str, str]] = {}
    for row in _read(genes_download.open(**kwargs), required=GENE_COLUMNS):
        concept = _cell(row['concept_id'])
        if concept is None:
            continue
        key = (row['gene_claim_name'], row['source_db_name'])
        value = (concept, _cell(row['gene_name']))
        if concepts.setdefault(key, value) != value:
            raise ValueError(f'DGIdb genes: the claim {key!r} has two concepts')
    return concepts


def iter_categories(
    opener,
    *,
    genes_download: Download,
    force_refresh: bool = False,
    **kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """Yield one record per gene category claim of ``categories.tsv``.

    A claim names a gene as its source does (a symbol, an alias, a UniProt or
    Ensembl accession, ``genes.tsv`` says which). The concept the claim maps to
    is joined from ``genes.tsv`` on the claim name and the source database.

    Output fields:

    - ``gene_claim_name``, ``category``, ``source``, ``source_version``: the
      columns ``gene_claim_name``, ``gene_category_name``, ``source_db_name``
      and ``source_db_version``.
    - ``gene_concept_id``: the concept of the claim in ``genes.tsv`` as written
      (``hgnc:1864``), None for a claim that ``genes.tsv`` maps to no gene.
    - ``gene_name``: the gene's name in ``genes.tsv`` (its HGNC symbol).
    """
    concepts = _gene_concepts(genes_download, force_refresh=force_refresh, **kwargs)
    for row in _read(opener, required=CATEGORY_COLUMNS):
        claim, source = row['gene_claim_name'], row['source_db_name']
        concept, name = concepts.get((claim, source), (None, None))
        yield {
            'gene_claim_name': claim,
            'category': _cell(row['gene_category_name']),
            'source': _cell(source),
            'source_version': _cell(row['source_db_version']),
            'gene_concept_id': concept,
            'gene_name': name,
        }


def iter_interactions(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield one record per interaction type of each drug-gene interaction claim.

    A claim without a type gives one record with ``interaction_type`` None. A
    claim with several types (``inhibitor|agonist``) gives one record per type,
    since each is a statement of the source. Rows of drugs that DGIdb did not
    map to a concept end after the source columns; their score columns are None.

    Output fields, the columns of the file with empty cells as None:

    - ``gene_claim_name``, ``gene_concept_id`` (``hgnc:10452``), ``gene_name``
      and ``drug_claim_name``, ``drug_concept_id`` (``chembl:CHEMBL592316``),
      ``drug_name``.
    - ``drug_is_approved``, ``drug_is_immunotherapy``, ``drug_is_antineoplastic``:
      ``true`` or ``false``.
    - ``source`` and ``source_version``: ``interaction_source_db_name`` and
      ``interaction_source_db_version``.
    - ``interaction_type``: one type, as written.
    - ``interaction_score``, ``drug_specificity_score``, ``gene_specificity_score``,
      ``evidence_score``: as written.
    """
    for row in _read(opener, required=INTERACTION_COLUMNS):
        record = {
            column: _cell(row[column])
            for column in INTERACTION_COLUMNS
            if column != 'interaction_types'
        }
        record['source'] = record.pop('interaction_source_db_name')
        record['source_version'] = record.pop('interaction_source_db_version')
        types = _cell(row['interaction_types'])
        for interaction_type in (types.split('|') if types else [None]):
            yield {**record, 'interaction_type': interaction_type}
