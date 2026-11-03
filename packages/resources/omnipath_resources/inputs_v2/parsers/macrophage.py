"""Read the interaction list of the macrophage pathway map.

Additional file 2 of Raza et al. 2010 lists the interactions drawn on the
authors' map of macrophage activation pathways. Each row names two interacting
partners by official gene symbol and Entrez Gene ID, with the kind of
interactant and its label on the map, the interaction type, its location and
one supporting PubMed ID; an interaction with several publications takes
several rows.

The worksheet ``Interaction Sheet Jan`` has two header rows, in which four
names occur once for each partner, so the columns are named by their partner
(``Partner 1 Official Gene Symbol``, ``Partner 2 Gene ID``, ...). The reader
yields one dictionary per interaction row with the twenty-one columns of
:data:`COLUMNS`: every value is text (a whole number without its decimal point,
as the sheet stores numbers as floats) or None for an empty cell. The
values of the ``Interactant Type`` columns must be among
:data:`INTERACTANT_TYPES`.
"""

from __future__ import annotations

from collections.abc import Generator, Iterable
import itertools
from typing import Any

from pypath.inputs_v2.parsers.base import iter_excel

# The worksheet with the published list; the workbook also has an older one.
SHEET = 'Interaction Sheet Jan'
# The two header rows of the sheet (rows 4 and 5). Four names occur once for
# each partner, so the shared worksheet reader cannot take the names from them.
HEADER_ROWS = (
    (
        '', 'Interacting Partner 1', '', '', '', 'Interacting Partner 2', '',
        '', '', 'Interaction Type', 'Interaction Location', 'NCBI-PubMed ID',
        'Species', 'Cell Type', 'Analytical Method', 'Additional Notes', '',
        '', '', '', '',
    ),
    (
        'Interaction No.', 'Official Gene Symbol', 'Gene ID',
        'Interactant Type', 'Interactant as on map', 'Official Gene Symbol',
        'Gene ID', 'Interactant Type', 'Interactant as on map', '', '', '',
        '', '', '', 'Interaction-description', 'Common aliases',
        'Other information', 'Pathway', 'Sub-Pathway', 'Contributor',
    ),
)
# The twenty-one columns of the sheet, in source order.
COLUMNS = (
    'Interaction No.',
    'Partner 1 Official Gene Symbol',
    'Partner 1 Gene ID',
    'Partner 1 Interactant Type',
    'Partner 1 Interactant as on map',
    'Partner 2 Official Gene Symbol',
    'Partner 2 Gene ID',
    'Partner 2 Interactant Type',
    'Partner 2 Interactant as on map',
    'Interaction Type',
    'Interaction Location',
    'NCBI-PubMed ID',
    'Species',
    'Cell Type',
    'Analytical Method',
    'Interaction-description',
    'Common aliases',
    'Other information',
    'Pathway',
    'Sub-Pathway',
    'Contributor',
)
# Every value of the two ``Interactant Type`` columns, as the sheet writes it.
# An empty cell is allowed. A value that is not here is a change of the sheet.
INTERACTANT_TYPES = frozenset(
    {
        'Protein',
        'Gene',
        'Complex',
        'Molecular Species',
        'dsRNA',
        'Drug',
        'Lipid:carbohydrate',
        'Generic Entity',
        'DNASeq [promotor element]',
        'DNASeq.',
        'N/A',
        'Phospholipid',
        'virus',
    }
)
_INTERACTANT_TYPE_COLUMNS = (
    'Partner 1 Interactant Type',
    'Partner 2 Interactant Type',
)


def _cell(value: Any) -> str | None:
    """A worksheet cell as text; a whole number without its decimal point."""
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value).strip() or None


def read_sheet(rows: Iterable[dict[str, Any]]) -> Generator[dict[str, Any], None, None]:
    """Yield the interactions of the sheet's rows, checking its two header rows.

    ``rows`` are the rows of the sheet from the fourth on, keyed by
    :data:`COLUMNS`: the two header rows, then the interactions. Empty cells
    are None.
    """
    rows = iter(rows)
    header = tuple(
        tuple('' if value is None else str(value) for value in row.values())
        for row in itertools.islice(rows, len(HEADER_ROWS))
    )
    if header != HEADER_ROWS:
        raise ValueError(
            f'Macrophage interaction sheet has unexpected header rows: {header}'
        )
    count = 0
    for row in rows:
        count += 1
        record = {
            column: None if value is None else _cell(value)
            for column, value in row.items()
        }
        # Checked on every row, whether or not the row is mapped.
        for column in _INTERACTANT_TYPE_COLUMNS:
            if record[column] not in INTERACTANT_TYPES | {None}:
                raise ValueError(
                    f'Unknown macrophage interactant type: {record[column]!r} '
                    f'({column}, interaction {record["Interaction No."]})'
                )
        yield record
    if not count:
        raise ValueError('Macrophage interaction sheet has no rows')


def iter_interactions(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield the rows of the interaction sheet of the workbook."""
    # The header cannot name the columns (see ``HEADER_ROWS``). Reading from
    # below row 3 gives its two rows first, so that they can be checked.
    yield from read_sheet(
        iter_excel(opener, sheet=SHEET, header_row=3, fieldnames=COLUMNS)
    )
