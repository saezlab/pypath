"""
Parser of the Guide to Pharmacology interaction table shipped with scConnect.

The file is `scConnect/data/GTP_tables_clean/2020-5/interactions.csv`, a comma
separated table with a header and 39 columns, one row for each measurement of a
ligand on a target. The parser yields every row keyed by the header; of the
columns these are read:

    target, target_gene_symbol, target_uniprot, target_ensembl_gene_id,
    target_species, ligand, ligand_gene_symbol, ligand_species, action, pubmed_id

A cell with several values joins them by `|`. `target_species` is one of
`Human`, `Mouse`, `Rat` (and others) for a target with a gene. A peptide ligand
lists the species it occurs in and its gene in each, but not in the same order
(`Human|Mouse|Rat` with `Avp|AVP|Avp`); `ligand_gene` reads the gene of the
target's species where the table allows it. A table without the columns, or
without rows, raises.
"""

from __future__ import annotations

from collections.abc import Generator
from typing import Any

from pypath.inputs_v2._publications import is_pmid
from pypath.inputs_v2.parsers.base import iter_csv

# The columns of `interactions.csv` that are read; the table has 39.
COLUMNS = (
    'target',
    'target_gene_symbol',
    'target_uniprot',
    'target_ensembl_gene_id',
    'target_species',
    'ligand',
    'ligand_gene_symbol',
    'ligand_species',
    'action',
    'pubmed_id',
)
# Every PubMed ID in the table has seven or eight digits.


def iter_interactions(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield the rows of the interaction table, checking that its columns are the expected ones."""
    count = 0
    for row in iter_csv(opener):
        missing = [column for column in COLUMNS if column not in row]
        if missing:
            raise ValueError(
                f'scConnect table lacks expected columns: {", ".join(missing)}'
            )
        count += 1
        yield row
    if not count:
        raise ValueError('scConnect table has no rows')


def split_values(value: Any) -> list[str]:
    """The values of a cell that separates them by `|`."""
    return [token for part in str(value or '').split('|') if (token := part.strip())]


def ligand_gene(row: dict[str, Any]) -> str | None:
    """The gene of a peptide ligand in the species of the row's target.

    The table lists the species a peptide occurs in and its gene in each, but
    not in the same order (`Human|Mouse|Rat` with `Avp|AVP|Avp`). The gene of
    the target's species is therefore known only when the peptide has one
    species and one gene, or the same symbol in each of its species.
    """
    genes = split_values(row.get('ligand_gene_symbol'))
    species = split_values(row.get('ligand_species'))
    if (
        len(genes) == len(species)
        and len(set(genes)) == 1
        and row.get('target_species') in species
    ):
        return genes[0]
    return None


def pubmed_ids(value: Any) -> list[str]:
    """The PubMed IDs of a `pubmed_id` cell, each once, in source order; a value
    that is not a PubMed ID is left out."""
    return [
        pmid
        for pmid in dict.fromkeys(split_values(value))
        if is_pmid(pmid)
    ]
