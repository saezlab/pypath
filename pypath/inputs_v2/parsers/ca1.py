"""Read the hippocampal CA1 neuron signaling network of Ma'ayan et al. 2005.

The network is file ``S1.txt`` of the supporting online material. Each line is
one link between two components with one supporting publication; a link with
several publications takes several lines. The columns are separated by single
spaces and four column names occur twice (``Human``, ``Mouse``, ``Function``
and ``Location``, once for the source and once for the target):

    Source Human Mouse Function Location Target Human Mouse Function Location Effects Type of interaction PubMedID
    CAMKII Q9UQM7 P11798 Kinase Cytosol SYNAPSIN P17600 O88935 Vesicle Vesicles _ Phosphorylation 7890745

The reader names the columns by their side (``Source Human``, ``Target
Function``, ...) and yields the thirteen columns of each line. ``NA`` (a
component without an accession: small molecules, ions, DNA elements and a few
proteins) is read as None. The ``Effects`` of every line must be one of
:data:`EFFECTS`. One line lacks a field; it is yielded with every
column None and the line itself as ``Malformed``.
"""

from __future__ import annotations

from collections.abc import Generator
import itertools
from typing import Any

from pypath.inputs_v2.parsers.base import iter_csv

# The header line of ``S1.txt``. The columns are separated by single spaces and
# four names occur twice, so the shared readers cannot take the names from it.
HEADER = (
    'Source Human Mouse Function Location '
    'Target Human Mouse Function Location '
    'Effects Type of interaction PubMedID'
)
# The thirteen columns of the table, in source order.
COLUMNS = (
    'Source',
    'Source Human',
    'Source Mouse',
    'Source Function',
    'Source Location',
    'Target',
    'Target Human',
    'Target Mouse',
    'Target Function',
    'Target Location',
    'Effects',
    'Type of interaction',
    'PubMedID',
)
# The table writes this for a component without an accession.
MISSING = 'NA'
# The values of the ``Effects`` column: activating, inhibitory and neutral.
EFFECTS = frozenset({'+', '_', '0'})


def _line(row: dict[Any, Any]) -> str:
    """The line of a row as the file writes it."""
    fields = [row[column] for column in COLUMNS if row[column] is not None]
    return ' '.join(fields + (row.get(None) or []))


def iter_links(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield the lines of ``S1.txt`` as rows, checking its header.

    ``NA`` becomes None. One line of the file lacks a field; it is yielded with
    every column None and the line itself as ``Malformed``.
    """
    # The header cannot name the columns (see ``HEADER``), so it is read as
    # the first row and compared as a line.
    rows = iter_csv(opener, delimiter=' ', fieldnames=COLUMNS, quotechar=None)
    header = [_line(row) for row in itertools.islice(rows, 1)]
    if not header or header[0].strip() != HEADER:
        raise ValueError(
            f'CA1 table does not start with the expected header: {header}'
        )
    count = 0
    for row in rows:
        line = _line(row)
        if not line.strip():
            continue
        count += 1
        if None in row or None in row.values():
            yield {**dict.fromkeys(COLUMNS), 'Malformed': line}
            continue
        # Checked on every line, whether or not the link is mapped.
        if row['Effects'] not in EFFECTS:
            raise ValueError(f'Unknown CA1 effect: {row["Effects"]!r}')
        yield {
            column: None if value == MISSING else value
            for column, value in row.items()
        }
    if not count:
        raise ValueError('CA1 table has no rows')
