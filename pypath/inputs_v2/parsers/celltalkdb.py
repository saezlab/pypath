"""
Read the CellTalkDB ligand-receptor tables.

Source files (R data frames, one per species):

    human_lr_pair.rds
    mouse_lr_pair.rds

A row is one ligand-receptor pair named by gene symbols, with NCBI Gene and
Ensembl cross-references for both partners and the PubMed IDs supporting the
pair (``evidence``, comma separated).
"""

from __future__ import annotations

from collections.abc import Generator
from typing import Any

from pypath.inputs_v2.parsers.base import iter_rdata

# Columns of ``human_lr_pair.rds`` and ``mouse_lr_pair.rds``, in source order.
COLUMNS = (
    'lr_pair',
    'ligand_gene_symbol',
    'receptor_gene_symbol',
    'ligand_gene_id',
    'receptor_gene_id',
    'ligand_ensembl_protein_id',
    'receptor_ensembl_protein_id',
    'ligand_ensembl_gene_id',
    'receptor_ensembl_gene_id',
    'evidence',
)


def iter_lr_pairs(
    opener,
    *,
    taxon_id: str,
    **_kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """
    Yield the rows of a CellTalkDB table with the file's species added.

    Every source column is kept under its own name. The table has no species
    column, so ``taxon_id`` is the NCBI Taxonomy ID of the file. Some cells
    hold the text ``NA`` for a missing cross-reference; they become None.
    """
    count = 0
    for row in iter_rdata(opener):
        missing = [column for column in COLUMNS if column not in row]
        if missing:
            raise ValueError(
                f'CellTalkDB table lacks expected columns: {", ".join(missing)}'
            )
        record = {
            column: None if value == 'NA' else value
            for column, value in row.items()
        }
        record['taxon_id'] = taxon_id
        count += 1
        yield record
    if not count:
        raise ValueError('CellTalkDB table has no rows')
