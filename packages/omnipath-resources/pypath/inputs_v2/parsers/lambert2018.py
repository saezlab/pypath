"""Read the table of human transcription factors of Lambert et al. 2018.

Table S1 of "The Human Transcription Factors" (Cell 172, 650) has one row per
candidate gene (2,765). The first row of its sheet groups the columns ("Gene
Information", "Is TF?", "Final Assesment", "Curator 1", "Prior
classifications" ...) and the second holds the column names, so the reader
names the columns by position and checks both header rows. The columns it
returns for each gene:

    id                Ensembl gene ID; ``uniprot_entry`` for the three rows
                      that have a UniProt entry name instead (see below)
    id_type           ``ensg`` or ``uniprot_entry``
    symbol            the gene's name in the table ("TF Name")
    dbd               the DNA-binding domain family that assigns the gene to
                      a TF family ("DBD"; ``Unknown`` where none is known)
    is_tf             the curated verdict, ``Yes`` or ``No``
    assessment        "TF assessment": the likelihood that the gene is a
                      sequence-specific TF, with emphasis on a known motif
    binding_mode      whether it binds alone or as a heteromer, as numbered
                      by the table (``1 Monomer or homomultimer`` ...)
    motif_status      whether a motif is known
    notes             "Notes" on the final annotation
    assessment_disagreement, binding_disagreement
                      ``Disagree`` where the two curators disagreed
    assessment_1, assessment_2, binding_1, binding_2
                      the assessment and binding mode of each curator

The source's ID column holds the UniProt entry name for three genes (DUX1,
DUX3 and ZNF73) where it should hold an Ensembl gene ID. The reader keeps
the value and names its kind in ``id_type``. The "Prior classifications"
columns (the classes of Vaquerizas 2009, CisBP, TFclass, TF-CAT and GO and
the PDB structures) are other resources' statements and are not read.
"""

from __future__ import annotations

import re
from collections.abc import Generator
from typing import Any

from pypath.inputs_v2.parsers.base import iter_excel

SHEET = 'Table S1. Related to Figure 1B'
# Names of the columns by position, as far as they are read.
FIELDS = (
    'id', 'symbol', 'dbd', 'is_tf', 'assessment', 'binding_mode',
    'motif_status', 'notes', 'comments', 'committee_notes', 'mtw_notes',
    'trh_notes', 'sl_notes', 'aj_notes', 'assessment_disagreement',
    'binding_disagreement', 'author_1', 'assessment_1', 'binding_1',
    'comment_1', 'notes_1', 'author_2', 'assessment_2', 'binding_2',
    'comment_2', 'notes_2', 'vaquerizas_2009', 'cisbp', 'tfclass',
    'tf_cat', 'go_tf', 'pdb', 'extra_1', 'extra_2', 'extra_3',
)
# What the second header row must say for the columns that are read, and what
# the first must say for the column that is named there.
SECOND_HEADER = {
    'id': 'ID',
    'symbol': 'Name',
    'dbd': 'DBD',
    'assessment': 'TF assessment',
    'binding_mode': 'Binding mode',
    'motif_status': 'Motif status',
    'notes': 'Notes',
    'assessment_disagreement': 'Disagree on Assessment',
    'binding_disagreement': 'Disagree on Binding',
    'assessment_1': 'Assesment1',
    'binding_1': 'Binding1',
    'assessment_2': 'Assesment2',
    'binding_2': 'Binding2',
}
FIRST_HEADER = {'is_tf': 'Is TF?'}
READ = (
    'id', 'symbol', 'dbd', 'is_tf', 'assessment', 'binding_mode',
    'motif_status', 'notes', 'assessment_disagreement',
    'binding_disagreement', 'assessment_1', 'assessment_2', 'binding_1',
    'binding_2',
)
UNKNOWN_DBD = 'Unknown'
# What joins the families of a gene with two kinds of DNA-binding domain.
DBD_SEPARATOR = '; '
_ENSEMBL_GENE = re.compile(r'^ENSG\d{11}$')
_UNIPROT_ENTRY = re.compile(r'^[A-Z0-9]+_HUMAN$')


def _check_header(first: dict[str, Any], second: dict[str, Any]) -> None:
    for header, expected in ((first, FIRST_HEADER), (second, SECOND_HEADER)):
        for field, name in expected.items():
            if (header.get(field) or '').strip() != name:
                raise ValueError(
                    f'Lambert 2018 Table S1: column {field!r} should be {name!r}, '
                    f'found {header.get(field)!r}'
                )


def iter_transcription_factors(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield one dictionary per gene of Table S1, with the fields above.

    Genes with the verdict ``No`` are included: they are statements that the
    gene is not a TF, which the mapping leaves out. A verdict that is neither
    ``Yes`` nor ``No`` and an ID that is neither an Ensembl gene ID nor a
    UniProt entry name are errors.
    """
    rows = iter_excel(opener, sheet=SHEET, header_row=0, fieldnames=FIELDS)
    first = next(rows, None)
    second = next(rows, None)
    if first is None or second is None:
        raise ValueError('Lambert 2018 Table S1 is empty')
    _check_header(first, second)
    count = 0
    for row in rows:
        gene = {field: row.get(field) for field in READ}
        if gene['is_tf'] not in {'Yes', 'No'}:
            raise ValueError(f'Unknown Lambert 2018 verdict in {gene!r}')
        if _ENSEMBL_GENE.match(str(gene['id'])):
            gene['id_type'] = 'ensg'
        elif _UNIPROT_ENTRY.match(str(gene['id'])):
            gene['id_type'] = 'uniprot_entry'
        else:
            raise ValueError(f'Malformed Lambert 2018 gene ID in {gene!r}')
        if not gene['symbol']:
            raise ValueError(f'Lambert 2018 gene without name: {gene!r}')
        count += 1
        yield gene
    if not count:
        raise ValueError('Lambert 2018 Table S1 has no genes')


def iter_dbd_families(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield one dictionary per DNA-binding domain family of the genes judged TFs.

    ``dbd`` is a family of the table and ``members`` the genes of the verdict
    ``Yes`` that have it, each with ``id``, ``id_type`` and ``symbol``. The table
    writes the families of a gene with two kinds of domain joined by ``; ``
    (``Homeodomain; POU``, 12 such pairs): the gene is a member of each of the
    two families, and the pair is no family of its own. Genes with the family
    ``Unknown`` are in no family. Families are in the order of first
    appearance.
    """
    families: dict[str, list[dict[str, Any]]] = {}
    for gene in iter_transcription_factors(opener):
        if gene['is_tf'] == 'Yes' and gene['dbd'] != UNKNOWN_DBD:
            if gene['dbd'] is None:
                raise ValueError(f'Lambert 2018 TF without DBD: {gene!r}')
            names = gene['dbd'].split(DBD_SEPARATOR)
            if any(not name.strip() or name != name.strip() for name in names) or len(
                set(names)
            ) != len(names):
                raise ValueError(f'Lambert 2018 DBD not understood: {gene["dbd"]!r}')
            for name in names:
                families.setdefault(name, []).append(
                    {key: gene[key] for key in ('id', 'id_type', 'symbol')}
                )
    for dbd, members in families.items():
        yield {'dbd': dbd, 'members': members}


def iter_list(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield one dictionary for the whole table: the number of genes by verdict.

    ``verdicts`` maps ``Yes`` and ``No`` to the number of genes. The table is
    read, so that a changed or empty file is an error.
    """
    verdicts = {'Yes': 0, 'No': 0}
    for gene in iter_transcription_factors(opener):
        verdicts[gene['is_tf']] += 1
    yield {'verdicts': verdicts}
