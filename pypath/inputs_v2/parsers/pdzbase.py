"""
Reader for the PDZBase "all interactions" page.

The tables of the page are read by the shared HTML reader;
:func:`iter_pdzbase` is the PDZBase interpretation of the page.
"""

from __future__ import annotations

from collections.abc import Generator
import re
from typing import Any

from pypath.inputs_v2.parsers.base import iter_html_tables, read_text

__all__ = ['iter_pdzbase']


# Column headers of the interactions table, in page order, with the field
# each is read into.
COLUMNS = (
    ('PDZ name', 'pdz_name'),
    ('PDZ accession code', 'pdz_accession'),
    ('PDZ domain', 'pdz_domain'),
    ('Ligand name', 'ligand_name'),
    ('Ligand accession code', 'ligand_accession'),
    ('Organism', 'organism'),
    ('Reference', 'reference'),
)
_HEADER = tuple(header for header, _ in COLUMNS)
# The last cell of a row links to the record's own page.
_DETAIL_LINK = re.compile(r'(?:^|/)interaction_detail/(\d+)$')


def iter_pdzbase(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield the interactions listed on the PDZBase "all interactions" page.

    The page holds one table whose first row is the seven column headers of
    :data:`COLUMNS`. Each following row is one interaction: the seven cells,
    a checkbox cell and a cell linking to ``interaction_detail/<number>``.

    Fields of a yielded record, all strings as written on the page:

    - ``pdz_name``, ``ligand_name``: PDZBase's names of the two proteins;
    - ``pdz_accession``, ``ligand_accession``: their UniProt accessions, some
      with an isoform number;
    - ``pdz_domain``: which PDZ domain of the protein binds, counted from the
      N-terminus;
    - ``organism``: PDZBase's organism label of the interaction;
    - ``reference``: the PubMed ID of the supporting publication;
    - ``interaction_id``: PDZBase's number of the record, from the detail
      link; None if the row has no such link.

    A page without the interactions table, a table without rows and a row
    with fewer than seven cells raise ``ValueError``.
    """
    html = read_text(opener)
    if not html.strip():
        return
    table = next(
        (
            table
            for table in iter_html_tables(html)
            if table['rows']
            and tuple(cell['text'] for cell in table['rows'][0]) == _HEADER
        ),
        None,
    )
    if table is None:
        raise ValueError(
            'PDZBase page has no table with the columns: ' + ', '.join(_HEADER)
        )
    rows = table['rows'][1:]
    if not rows:
        raise ValueError('PDZBase interactions table has no rows')
    for number, cells in enumerate(rows, 1):
        if len(cells) < len(COLUMNS):
            raise ValueError(
                f'PDZBase interaction row {number} has {len(cells)} cells, '
                f'expected at least {len(COLUMNS)}'
            )
        record = {
            field: cell['text']
            for (_, field), cell in zip(
                COLUMNS, cells[: len(COLUMNS)], strict=True
            )
        }
        record['interaction_id'] = next(
            (
                match.group(1)
                for cell in cells[len(COLUMNS):]
                for link in cell['links']
                if (match := _DETAIL_LINK.search(link))
            ),
            None,
        )
        yield record
