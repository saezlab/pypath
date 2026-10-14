"""
Read the miRDeathDB table of microRNAs and their targets in programmed cell death.

Source file (tab separated, with a header row, fields with commas in double quotes):

    miRDeathDB_all_data.txt

A row is one curated observation: a microRNA (``miRNA_symbol``) and a target gene
(``Gene_Symbol``, Entrez ``geneid``) in a cell death pathway (``Pathway``), the
organism and tissue, the PubMed ID and the nomenclature of the gene. The
``miRBase_mature_ID`` and ``miRBase_ID`` cells list, comma separated and in the same
order, the miRBase mature names and accessions (MIMAT) that the curators matched to
the symbol; a symbol such as ``miR-221`` lists both arms of the hairpin, one such as
``miR-106`` lists several hairpins.

``Action_Mode`` is the expression direction of the microRNA in the experiment (the
ncRDeathDB detail page of the same data names it "Expression Direction of ncRNA").

The helpers at the top (``split_cell``, ``mirbase_products``, ``gene_links``) also
serve the ncRDeathDB parser, whose rows are written the same way.
"""

from __future__ import annotations

import re
from collections.abc import Generator
from typing import Any

from pypath.inputs_v2._publications import is_pmid
from pypath.inputs_v2.parsers.base import iter_tsv

# Columns of ``miRDeathDB_all_data.txt``, in source order.
COLUMNS = (
    'miRNA_symbol',
    'miRBase_mature_ID',
    'miRBase_ID',
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
)

# Cells that stand for a missing value.
_NULLS = frozenset({'', '-', 'NA', 'N/A'})
_mimat = re.compile(r'MIMAT\d{7}')
# A miRBase mature name: the organism code, then the name (``hsa-miR-21-5p``).
_mature_name = re.compile(r'([a-z]{3,4})-(?:miR|mir|let|lsy|bantam)[A-Za-z0-9*.\-]*')
_entrez = re.compile(r'\d+')
_ensg = re.compile(r'ENSG\d{11}')


def cell(value: Any) -> str | None:
    """A cell stripped of white space; None for the source's missing-value marks."""
    value = (value or '').strip()
    return None if value in _NULLS else value


def split_cell(value: Any, separator: str = '|') -> list[str]:
    """The non-empty parts of a list cell, stripped, in source order."""
    return [part.strip() for part in (value or '').split(separator) if cell(part)]


def mirbase_products(names: Any, accessions: Any) -> dict[str, Any]:
    """The miRBase products the curators matched to a microRNA symbol.

    ``names`` and ``accessions`` are the comma separated cells. ``mirbase_mature``
    and ``mirbase_mature_name`` are set only when the symbol has exactly one
    well-formed product: one MIMAT accession and one mature name. A list of
    products is kept as it is, with the organism codes of its names.
    """
    name_list = split_cell(names, ',')
    accession_list = split_cell(accessions, ',')
    single = (
        len(name_list) == 1
        and len(accession_list) == 1
        and _mimat.fullmatch(accession_list[0])
        and _mature_name.fullmatch(name_list[0])
    )
    return {
        'mature_names': name_list,
        'mature_accessions': accession_list,
        'mirbase_mature': accession_list[0] if single else None,
        'mirbase_mature_name': name_list[0] if single else None,
        'mature_organism_codes': sorted({
            match.group(1) for name in name_list
            if (match := _mature_name.match(name))
        }),
    }


def each_product(products: dict[str, Any]) -> list[dict[str, Any]]:
    """The fields of ``mirbase_products`` once for each MIMAT accession listed.

    As legacy pypath, a symbol with several products gives an observation for
    each of them, and so does a cell with a malformed accession beside a valid
    one (``MIMAT0000257,MIMAT00226921``; the malformed one is left out). The
    mature name is that of the same position where the cell of names lists as
    many names as accessions.
    """
    accessions = [a for a in products['mature_accessions'] if _mimat.fullmatch(a)]
    if products['mirbase_mature'] or not accessions:
        return [products]
    names = products['mature_names']
    aligned = len(names) == len(products['mature_accessions'])
    result = []
    for accession in accessions:
        name = names[products['mature_accessions'].index(accession)] if aligned else None
        name = name if name and _mature_name.fullmatch(name) else None
        result.append({
            **products,
            'mirbase_mature': accession,
            'mirbase_mature_name': name,
            'mature_organism_codes': (
                [_mature_name.match(name).group(1)] if name else products['mature_organism_codes']
            ),
        })
    return result


def gene_links(value: Any) -> dict[str, Any]:
    """The HGNC ID and the Ensembl human gene ID of a ``Links`` cell, where it has exactly one."""
    links = {}
    for link in split_cell(value):
        database, _, accession = link.partition(':')
        links.setdefault(database, []).append(accession.strip())
    hgnc = [a for a in links.get('HGNC', []) if a.isdigit()]
    ensg = [a for a in links.get('Ensembl', []) if _ensg.fullmatch(a)]
    return {
        'hgnc': hgnc[0] if len(hgnc) == 1 else None,
        'ensembl_gene': ensg[0] if len(ensg) == 1 else None,
        'link_databases': sorted(links),
    }


def gene_fields(row: dict[str, Any]) -> dict[str, Any]:
    """The target gene of a miRDeathDB or ncRDeathDB row."""
    entrez = cell(row.get('geneid'))
    if entrez is not None and not _entrez.fullmatch(entrez):
        raise ValueError(f'Entrez gene ID is not a number: {entrez!r}')
    return {
        'gene_symbol': cell(row.get('Gene_Symbol')),
        'entrez_id': entrez,
        'gene_synonyms': split_cell(cell(row.get('Synonyms'))),
        'gene_name': cell(row.get('Full_name_from_nomenclature_authority'))
        or cell(row.get('Description')),
        'gene_type': cell(row.get('type_of_gene')),
        **gene_links(row.get('Links')),
    }


def iter_interactions(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield the observations of miRDeathDB, checking that the columns are the expected ones.

    Fields: ``rna_symbol`` (the microRNA symbol of the source), ``mature_names``,
    ``mature_accessions``, ``mirbase_mature`` and ``mirbase_mature_name`` (one
    observation for each product of the symbol, see ``each_product``),
    ``mature_organism_codes``; the gene fields of ``gene_fields``; ``pathway``,
    ``expression_direction`` (``Action_Mode``), ``organism``, ``tax_id`` (as in the
    source), ``tissue``; ``pmids`` (list of PubMed IDs; the PMID cell as a
    list when it is a PubMed ID).
    """
    count = 0
    rows = iter_tsv(opener, quotechar='"', fieldnames=COLUMNS, strict=True)
    header = next(rows, None)
    if header is None:
        raise ValueError('miRDeathDB table is empty')
    if tuple(header.values()) != COLUMNS:
        raise ValueError(
            f'miRDeathDB header differs from the expected columns: {tuple(header.values())}'
        )
    for row in rows:
        count += 1
        pmid = cell(row['PMID'])
        if pmid is not None and not pmid.isdigit():
            raise ValueError(f'miRDeathDB PMID is not a number: {pmid!r}')
        for products in each_product(mirbase_products(row['miRBase_mature_ID'], row['miRBase_ID'])):
            yield {
                'rna_symbol': cell(row['miRNA_symbol']),
                **products,
                **gene_fields(row),
                'pathway': cell(row['Pathway']),
                'expression_direction': cell(row['Action_Mode']),
                'organism': cell(row['Organism']),
                'tax_id': cell(row['tax_id']),
                'tissue': cell(row['Tissue']),
                'pmids': [pmid] if is_pmid(pmid) else [],
            }
    if not count:
        raise ValueError('miRDeathDB table has no rows')
