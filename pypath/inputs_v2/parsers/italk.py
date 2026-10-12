"""
Read the iTALK ligand-receptor database.

Source file (R data file of the iTALK package, object ``database``):

    LR_database.rda

A row is one human ligand-receptor pair: ``Pair.Name``, the approved gene
symbol and the name of the ligand and of the receptor, and a
``Classification`` (cytokine, growth factor, checkpoint or other). A symbol
cell is empty for the viral receptors, which the source names only by alias.
"""

from __future__ import annotations

from collections.abc import Generator
from typing import Any

from pypath.inputs_v2.parsers.base import iter_rdata

# Columns of the ``database`` data frame, in source order.
COLUMNS = (
    'Pair.Name',
    'Ligand.ApprovedSymbol',
    'Ligand.Name',
    'Receptor.ApprovedSymbol',
    'Receptor.Name',
    'Classification',
)


def iter_pairs(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield the rows of iTALK's data frame, checking that its columns are the expected ones."""
    count = 0
    for row in iter_rdata(opener, name='database'):
        missing = [column for column in COLUMNS if column not in row]
        if missing:
            raise ValueError(
                f'iTALK table lacks expected columns: {", ".join(missing)}'
            )
        count += 1
        yield row
    if not count:
        raise ValueError('iTALK table has no rows')
