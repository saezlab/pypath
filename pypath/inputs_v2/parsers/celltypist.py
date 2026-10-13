"""
Parser of the CellTypist encyclopedia table of immune cell types.

``encyclopedia_table.xlsx`` (worksheet ``Main``) has one row for each cell type
of the Pan-immune CellTypist atlas, version 2. Its columns are:

    High-hierarchy cell types       the broader group, e.g. ``B cells``
    Low-hierarchy cell types        the cell type itself, e.g. ``Follicular B cells``
    Description                     one sentence about the cell type
    Cell Ontology ID                the Cell Ontology term of the cell type
                                    (``CL:0000843``); several cell types of the
                                    table can have the same term
    Tissues                         the tissues the cell type was found in, ``, `` separated
    Datasets                        the datasets it was found in, ``, `` separated
    Curated markers                 marker genes chosen by the curators, ``, `` separated
    Top 10 important genes from the CellTypist model
                                    the ten genes with the largest weight for the
                                    cell type in the CellTypist model, ``, `` separated

:func:`iter_markers` yields one row for each gene of the two gene columns.
"""

from __future__ import annotations

from collections.abc import Generator
import functools
import re
from typing import Any

from pypath.inputs_v2.parsers.base import iter_excel

CURATED_MARKERS = 'Curated markers'
MODEL_GENES = 'Top 10 important genes from the CellTypist model'
COLUMNS = (
    'High-hierarchy cell types',
    'Low-hierarchy cell types',
    'Description',
    'Cell Ontology ID',
    'Tissues',
    'Datasets',
    CURATED_MARKERS,
    MODEL_GENES,
)
# The column that a gene comes from, as the ``marker_source`` of a row.
GENE_COLUMNS = {
    'curated_marker': CURATED_MARKERS,
    'model_top_gene': MODEL_GENES,
}

_CL_ID = re.compile(r'CL:\d{7}')
_LIST_SEPARATOR = re.compile(r'\s*,\s*')

_read = functools.partial(iter_excel, sheet='Main', header_row=1, required=COLUMNS)


def _items(cell: str | None) -> list[str]:
    """The items of a ``, `` separated cell, in order and without repeats."""
    if not cell:
        return []
    return list(dict.fromkeys(item for item in _LIST_SEPARATOR.split(cell.strip()) if item))


def iter_markers(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """One row for each gene of the curated markers and of the model genes.

    Fields: ``gene`` (symbol as written), ``marker_source``
    (``curated_marker`` or ``model_top_gene``), ``cell_type`` (low
    hierarchy), ``cell_type_group`` (high hierarchy), ``cell_ontology_id``,
    ``description``, ``tissues`` and ``datasets`` (lists). A cell type without a
    well formed Cell Ontology ID raises: the table gives every cell type one.
    """
    for row in _read(opener):
        cell_ontology_id = row['Cell Ontology ID']
        if not cell_ontology_id or not _CL_ID.fullmatch(cell_ontology_id):
            raise ValueError(
                f'CellTypist cell type {row["Low-hierarchy cell types"]!r} '
                f'has no Cell Ontology ID: {cell_ontology_id!r}'
            )
        common = {
            'cell_type': row['Low-hierarchy cell types'],
            'cell_type_group': row['High-hierarchy cell types'],
            'cell_ontology_id': cell_ontology_id,
            'description': row['Description'],
            'tissues': _items(row['Tissues']),
            'datasets': _items(row['Datasets']),
        }
        for source, column in GENE_COLUMNS.items():
            for gene in _items(row[column]):
                yield {**common, 'gene': gene, 'marker_source': source}
