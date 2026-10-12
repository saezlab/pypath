"""
Read the talklr ligand-receptor database.

Source file (R data file of the talklr package, object ``receptor_ligand``):

    receptor_ligand.rda

A row is one human ligand-receptor pair of Ramilowski et al. 2015: the pair
name, the approved gene symbol and name of ligand and receptor, one column per
database that found the pair (the cell repeats the column name or is empty),
the PubMed IDs the authors added (``PMID.Manual``, comma separated), and the
verdict of the curation in ``Pair.Evidence``: ``literature supported`` or
``putative``. The reader selects the rows of one verdict; a database cell that is
neither empty nor the column's name raises.
"""

from __future__ import annotations

from collections.abc import Generator
from typing import Any

from pypath.inputs_v2.parsers.base import iter_rdata

# Columns of the ``receptor_ligand`` data frame, in source order.
COLUMNS = (
    'Pair.Name',
    'Ligand.ApprovedSymbol',
    'Ligand.Name',
    'Receptor.ApprovedSymbol',
    'Receptor.Name',
    'DLRP',
    'HPMR',
    'IUPHAR',
    'HPRD',
    'STRING.binding',
    'STRING.experiment',
    'HPMR.Ligand',
    'HPMR.Receptor',
    'PMID.Manual',
    'Pair.Source',
    'Pair.Evidence',
)
# The database columns: a cell repeats the column's name or is empty.
DATABASE_COLUMNS = (
    'DLRP',
    'HPMR',
    'IUPHAR',
    'HPRD',
    'STRING.binding',
    'STRING.experiment',
)
# Every value of ``Pair.Evidence``, the verdict of the curation: a pair with
# primary literature support, or one found only in a database.
EVIDENCE = ('literature supported', 'putative')


def iter_pairs(
    opener,
    *,
    evidence: str,
    **_kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """
    Yield the rows of talklr's data frame that carry one ``Pair.Evidence`` value.

    Every source column is kept under its own name. The columns are checked
    to be the expected ones, the verdict of every row to be a known one and
    every database flag to be the column's name or empty.
    """
    count = 0
    for row in iter_rdata(opener, name='receptor_ligand'):
        missing = [column for column in COLUMNS if column not in row]
        if missing:
            raise ValueError(
                f'talklr table lacks expected columns: {", ".join(missing)}'
            )
        if row['Pair.Evidence'] not in EVIDENCE:
            raise ValueError(
                f'Unknown talklr pair evidence {row["Pair.Evidence"]!r}'
            )
        for column in DATABASE_COLUMNS:
            if row[column] not in (None, column):
                raise ValueError(
                    f'Unknown talklr database flag {row[column]!r} in {column}'
                )
        count += 1
        if row['Pair.Evidence'] == evidence:
            yield row
    if not count:
        raise ValueError('talklr table has no rows')
