"""
Read the LRdb ligand-receptor table.

Source file (tab separated, with a header row):

    LRdb_122019.txt

A row is one human ligand-receptor pair: the gene symbols (``ligand``,
``receptor``), protein names and synonyms for both, the databases the pair
was taken from (``source``, comma separated), the PubMed IDs (``PMIDs``,
comma separated), the cell types of ligand and receptor, and a free-text
``remarks`` cell. The rows are yielded as the file has them. A value of
``source`` that is not one of the known ones raises.
"""

from __future__ import annotations

from collections.abc import Generator
from typing import Any

from pypath.inputs_v2.parsers.base import iter_tsv

# Columns of ``LRdb_122019.txt``, in source order.
COLUMNS = (
    'ligand',
    'receptor',
    'ligand.name',
    'receptor.name',
    'ligand.synonyms',
    'receptor.synonyms',
    'ligand.altern.names',
    'receptor.altern.names',
    'source',
    'PMIDs',
    'cells.L',
    'cells.R',
    'remarks',
)
# Every value of the ``source`` column: seven databases, in the source's
# spelling, and ``literature``, a pair the LRdb authors curated themselves.
SOURCES = (
    'fantom5',
    'HPRD',
    'HPMR',
    'IUPHAR',
    'reactome',
    'uniprot',
    'cellsignal.com',
    'literature',
)


def iter_pairs(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield the rows of the LRdb table.

    The columns are checked to be the expected ones and every comma separated
    value of ``source`` to be a known one.
    """
    count = 0
    for row in iter_tsv(opener):
        missing = [column for column in COLUMNS if column not in row]
        if missing:
            raise ValueError(
                f'LRdb table lacks expected columns: {", ".join(missing)}'
            )
        for value in row['source'].split(','):
            if value not in SOURCES:
                raise ValueError(
                    f'Unknown LRdb source {value!r} for {row["ligand"]} - '
                    f'{row["receptor"]}'
                )
        count += 1
        yield row
    if not count:
        raise ValueError('LRdb table has no rows')
