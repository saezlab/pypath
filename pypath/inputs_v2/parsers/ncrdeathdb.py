"""
Read the ncRDeathDB table of non-coding RNAs in programmed cell death.

Source file (tab separated, ISO-8859-1, with a header row, no quoting):

    allNcRNACelldeathData.txt

A row is one observation: a non-coding RNA (``miRNA_symbol``) of a category
(``RNA Category``: miRNA, lncRNA or snoRNA), an interacting gene (``Gene_Symbol``,
Entrez ``geneid``; empty when the entry only links the RNA to a death pathway), the
death pathway, the organism and tissue, the PubMed ID and a description of the
evidence. The gene columns are those of miRDeathDB. ``Action_Mode`` is the expression
direction of the RNA in the experiment ("Expression Direction of ncRNA" on the
database's detail pages).

For a microRNA, ``miRBase_mature_ID`` and ``miRBase_ID`` list the miRBase products
that the curators matched to the symbol; several rows write a family instead
("miR-17/17-5p/20ab/...") or give lists that differ in length or order. For a
lncRNA or snoRNA the ``miRBase_mature_ID`` cell holds a link to a web page.

Rows whose PMID cell reads ``prediction`` are homology predictions: the entry of
another species is predicted from a curated human entry (the description names it).
They are separate from the curated observations, and each selection has its own
reader.
"""

from __future__ import annotations

from collections.abc import Callable, Generator
from typing import Any

from pypath.inputs_v2._publications import is_pmid
from pypath.inputs_v2.parsers.base import iter_tsv
from pypath.inputs_v2.parsers.mirdeathdb import (
    cell,
    gene_fields,
    mirbase_products,
)

# Column names of ``allNcRNACelldeathData.txt``, in source order. The header
# writes ``Description`` twice, for the gene and for the evidence; the second is
# called ``Evidence`` here.
SOURCE_HEADER = (
    'id',
    'miRNA_symbol',
    'miRBase_mature_ID',
    'miRBase_ID',
    'RNA Category',
    'Gene_Symbol',
    'Pathway',
    'Action_Mode',
    'Organism',
    'Tissue',
    'PMID',
    'tax_id',
    'geneid',
    'Synonyms',
    'Links',
    'chromosome',
    'map_location',
    'Description',
    'type_of_gene',
    'Full_name_from_nomenclature_authority',
    'Other_designations',
    'Description',
)
COLUMNS = SOURCE_HEADER[:-1] + ('Evidence',)

# The PMID cell of a homology prediction.
PREDICTION = 'prediction'


def _observations(opener) -> Generator[dict[str, Any], None, None]:
    """Every row of the table with its fields split and cleaned."""
    rows = iter_tsv(opener, quotechar=None, fieldnames=COLUMNS, strict=True)
    header = next(rows, None)
    if header is None:
        raise ValueError('ncRDeathDB table is empty')
    if tuple(header.values()) != SOURCE_HEADER:
        raise ValueError(
            f'ncRDeathDB header differs from the expected columns: {tuple(header.values())}'
        )
    for row in rows:
        pmid = cell(row['PMID'])
        if pmid is not None and pmid != PREDICTION and not pmid.isdigit():
            raise ValueError(f'ncRDeathDB PMID is neither a number nor "prediction": {pmid!r}')
        category = cell(row['RNA Category'])
        yield {
            'id': cell(row['id']),
            'rna_symbol': cell(row['miRNA_symbol']),
            'rna_category': category,
            **(
                mirbase_products(row['miRBase_mature_ID'], row['miRBase_ID'])
                if category == 'miRNA'
                else {}
            ),
            **gene_fields(row),
            'pathway': cell(row['Pathway']),
            'expression_direction': cell(row['Action_Mode']),
            'organism': cell(row['Organism']),
            'tax_id': cell(row['tax_id']),
            'tissue': cell(row['Tissue']),
            'pmids': [pmid] if is_pmid(pmid) else [],
            'prediction': pmid == PREDICTION,
            'description': cell(row['Evidence']),
        }


def _selection(prediction: bool, label: str) -> Callable[..., Generator[dict[str, Any], None, None]]:
    def iterate(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
        count = 0
        for row in _observations(opener):
            if row['prediction'] is prediction:
                count += 1
                yield row
        if not count:
            raise ValueError(f'ncRDeathDB table has no {label}')

    iterate.__doc__ = (
        f'Yield the {label} of ncRDeathDB, checking that the columns are the expected ones.\n\n'
        'Fields: ``id``, ``rna_symbol``, ``rna_category`` (miRNA, lncRNA or snoRNA, '
        'stripped), for a miRNA the miRBase fields of ``mirbase_products``, the gene '
        'fields of ``gene_fields``, ``pathway``, ``expression_direction``, ``organism``, '
        '``tax_id``, ``tissue``, ``pmids`` (PubMed IDs), ``prediction`` and '
        '``description`` (the prose of the evidence).'
    )
    return iterate


iter_interactions = _selection(False, 'curated observations')
iter_homology_predictions = _selection(True, 'homology predictions')
