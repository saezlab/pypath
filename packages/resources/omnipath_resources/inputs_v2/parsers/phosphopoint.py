"""Read the PhosphoPOINT human kinase interactome table.

The table pairs each kinase with a protein and puts the pair into one of four
categories (semicolon separated, with a header):

    Kinase;Entrez GeneID;Gene;Entrez GeneID;Category
    AKT1;207;ESR2;2100;Category 2
"""

from __future__ import annotations

from collections.abc import Generator
from typing import Any

from pypath.inputs_v2.parsers.base import iter_semicolon

# The header line. ``Entrez GeneID`` occurs twice, so the shared readers
# cannot take the column names from it.
HEADER = ('Kinase', 'Entrez GeneID', 'Gene', 'Entrez GeneID', 'Category')
# The five columns of the table, in source order.
COLUMNS = (
    'Kinase',
    'Kinase Entrez GeneID',
    'Gene',
    'Gene Entrez GeneID',
    'Category',
)


def iter_pairs(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield the rows of the PhosphoPOINT table, checking its header."""
    # The header cannot name the columns (see ``HEADER``), so it is read as
    # the first row. A row with another number of fields is an error.
    rows = (
        {column: value.strip() for column, value in row.items()}
        for row in iter_semicolon(opener, fieldnames=COLUMNS, strict=True)
    )
    header = tuple(next(rows, {}).values())
    if header != HEADER:
        raise ValueError(f'PhosphoPOINT table has an unexpected header: {header}')
    count = 0
    for row in rows:
        count += 1
        yield row
    if not count:
        raise ValueError('PhosphoPOINT table has no rows')
