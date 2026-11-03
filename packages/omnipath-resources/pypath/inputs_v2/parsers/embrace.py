"""
Read the EMBRACE ligand-receptor table.

Source file (Excel workbook, Table S11 of Sheikh et al. 2019, sheet
``expressed_pairs``):

    1-s2.0-S2589004219304055-mmc2.xlsx

A row is one mouse ligand-receptor pair named by gene symbols
(``ligand_symbol``, ``receptor_symbol``), followed by eight flags saying in
which of four cell types of the embryonic brain the ligand and the receptor
are expressed.
"""

from __future__ import annotations

from collections.abc import Generator
from typing import Any

from pypath.inputs_v2.parsers.base import iter_excel

SHEET = 'expressed_pairs'
# Columns of the sheet, in source order. The eight flags are 1 when the
# ligand or the receptor is expressed in neural (N), mural (P), microglial
# (M) or endothelial (EC) cells and 0 otherwise.
COLUMNS = (
    'ligand_symbol',
    'receptor_symbol',
    'N_ligand',
    'P_ligand',
    'M_ligand',
    'EC_ligand',
    'N_receptor',
    'P_receptor',
    'M_receptor',
    'EC_receptor',
)


def iter_pairs(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield the rows of Table S11, checking that its columns are the expected ones."""
    count = 0
    for row in iter_excel(opener, sheet=SHEET):
        missing = [column for column in COLUMNS if column not in row]
        if missing:
            raise ValueError(
                f'EMBRACE table lacks expected columns: {", ".join(missing)}'
            )
        count += 1
        yield row
    if not count:
        raise ValueError('EMBRACE table has no rows')
