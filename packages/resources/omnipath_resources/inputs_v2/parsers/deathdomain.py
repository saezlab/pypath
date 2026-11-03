"""
Parser of the DeathDomain database pages.

The database had one HTML page per subfamily (CARD, DD, DED, PYD), each with
four tables: "In vitro interaction", "In vivo interaction", "Characterization"
and "Functional Role". The tables are read with the shared table reader
(:func:`iter_html_tables`). A table line is one observation from one
publication: two proteins, a method, the details of that table and a PubMed ID.

The database lists every observation between two different proteins twice,
once from each protein; :func:`mark_second_listings` marks the second listing.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Generator
from typing import Any

from pypath.inputs_v2.parsers.base import iter_html_tables, read_text

# Table title -> (table name, columns of a data line). The third cell of a
# line holds only a link to the page of the pair.
_TABLES = {
    'In vitro interaction': (
        'in_vitro',
        (
            'family', 'protein_a', 'pair_link', 'protein_b', 'method',
            'a_species', 'a_region', 'a_expression',
            'b_species', 'b_region', 'b_expression',
            'reference',
        ),
    ),
    'In vivo interaction': (
        'in_vivo',
        (
            'family', 'protein_a', 'pair_link', 'protein_b', 'method',
            'endogenous_expression', 'overexpression',
            'a_species', 'a_region', 'b_species', 'b_region',
            'reference',
        ),
    ),
    'Characterization': (
        'characterization',
        (
            'family', 'protein_a', 'pair_link', 'protein_b', 'method',
            'stoichiometry', 'affinity', 'mutation', 'complex_structure',
            'reference',
        ),
    ),
    'Functional Role': (
        'functional_role',
        (
            'family', 'protein_a', 'pair_link', 'protein_b', 'method',
            'death_related', 'death_unrelated',
            'reference',
        ),
    ),
}
# The last header line of each table, which names the leading columns.
_HEADER_START = ['Family', 'DD1', '', 'DD2', 'Method']


def iter_family_page(html: str, family: str) -> Generator[dict[str, str], None, None]:
    """Yield the data lines of the four tables of one subfamily page.

    Each line is a dictionary of the table's columns (see ``_TABLES``) with
    ``table``, the name of the table, and ``page_family``. A page without the
    four tables, or a table with another layout, raises ``ValueError``.
    """
    found = set()
    for table in iter_html_tables(html):
        if table['attrs'].get('class') != 'tab' or not table['rows']:
            continue
        title = table['rows'][0][0]['text']
        if title not in _TABLES:
            raise ValueError(f'DeathDomain {family}: unknown table {title!r}')
        name, columns = _TABLES[title]
        found.add(name)
        headers = [
            [cell['text'] for cell in row]
            for row in table['rows']
            if all(cell['header'] for cell in row)
        ]
        if headers[-1][: len(_HEADER_START)] != _HEADER_START:
            raise ValueError(
                f'DeathDomain {family}: unexpected header of {title!r}: '
                f'{headers[-1]}'
            )
        for row in table['rows']:
            if any(cell['header'] for cell in row):
                continue
            if len(row) != len(columns):
                raise ValueError(
                    f'DeathDomain {family}: a line of {title!r} has '
                    f'{len(row)} cells, expected {len(columns)}'
                )
            record = {column: cell['text'] for column, cell in zip(columns, row)}
            del record['pair_link']
            yield {'table': name, 'page_family': family, **record}
    missing = {name for name, _ in _TABLES.values()} - found
    if missing:
        raise ValueError(
            f'DeathDomain {family}: page lacks the tables {sorted(missing)}'
        )


def _mirror(row: dict[str, str]) -> tuple[tuple[str, str], ...]:
    """The line as the page lists it from the other protein."""
    swapped = dict(row)
    for column in row:
        if column.startswith('a_') or column == 'protein_a':
            other = column.replace('a_', 'b_', 1) if column != 'protein_a' else 'protein_b'
            swapped[column], swapped[other] = row[other], row[column]
    return tuple(sorted(swapped.items()))


def mark_second_listings(rows) -> Generator[dict[str, Any], None, None]:
    """Add ``mirrored``: True for the second of two lines that differ only by
    which protein is listed first."""
    awaited: Counter = Counter()
    for row in rows:
        key = tuple(sorted(row.items()))
        mirrored = False
        if row['protein_a'] != row['protein_b']:
            if awaited[key]:
                awaited[key] -= 1
                mirrored = True
            else:
                awaited[_mirror(row)] += 1
        yield {**row, 'mirrored': mirrored}


def iter_deathdomain_observations(
    opener=None,
    force_refresh: bool = False,
    family_downloads: dict[str, Any] | None = None,
    **_kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """Yield one dictionary per table line of the four subfamily pages.

    ``family_downloads`` maps each subfamily to the download of its page.

    Output fields:

    - ``table``: ``in_vitro``, ``in_vivo``, ``characterization`` or
      ``functional_role``.
    - ``page_family``: the subfamily page; ``family``: the line's own cell.
    - ``protein_a``, ``protein_b``: protein names as the database writes them.
    - ``method``: the method label.
    - ``reference``: PubMed ID.
    - the other columns of the table, as written: ``a_species``, ``a_region``,
      ``a_expression`` and the same for ``b`` (in vitro); ``a_species``,
      ``a_region``, ``b_species``, ``b_region``, ``endogenous_expression``,
      ``overexpression`` (in vivo); ``stoichiometry``, ``affinity``,
      ``mutation``, ``complex_structure`` (characterization);
      ``death_related``, ``death_unrelated`` (functional role).
    - ``mirrored``: True for the second listing of an observation.
    """
    if not family_downloads:
        raise ValueError('DeathDomain: the downloads of the pages are needed')
    for family, download in family_downloads.items():
        html = read_text(download.open(force_refresh=force_refresh))
        yield from mark_second_listings(iter_family_page(html, family))
