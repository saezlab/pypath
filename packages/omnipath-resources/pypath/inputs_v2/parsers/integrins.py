"""Read Table 1 of "The integrins" (Takada et al. 2007, Genome Biology).

The table lists the human integrin subunits:

    Gene symbol  Protein name  Synonyms  Gene accession number
    Location (chromosome)  Protein accession number
    ITGA1        α1            CD49a     NM_181501              5q11.2  P56199

The article is read from Europe PMC as JATS XML: PubMed Central now asks a
client that is not a browser to solve a CAPTCHA before it serves the article's
table pages.
"""

from __future__ import annotations

from collections.abc import Generator
from typing import Any
import xml.etree.ElementTree as ET

from pypath.inputs_v2.parsers.base import read_xml_root

TABLE_ID = 'T1'
CAPTION = 'Human integrin subunits'
COLUMNS = (
    'Gene symbol',
    'Protein name',
    'Synonyms',
    'Gene accession number',
    'Location (chromosome)',
    'Protein accession number',
)


def _text(element: ET.Element) -> str:
    return ' '.join(''.join(element.itertext()).split())


def iter_subunits(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """One record per subunit of the table, in the table's order: ``symbol``,
    ``name`` (the subunit's name, ``α1``), ``synonyms`` (the cell, None if
    empty), ``gene_accession``, ``location`` and ``accession`` (of the protein:
    a UniProt accession, or one of GenBank where the table has no UniProt one)."""
    root = read_xml_root(opener)
    if root is None:
        raise ValueError('No integrin article was downloaded')
    wrap = next(
        (w for w in root.iter('table-wrap') if w.get('id') == TABLE_ID), None
    )
    if wrap is None:
        raise ValueError(f'The article has no table {TABLE_ID}')
    caption = _text(wrap.find('caption'))
    if caption != CAPTION:
        raise ValueError(f'Table {TABLE_ID} is {caption!r}, not {CAPTION!r}')
    rows = [[_text(cell) for cell in tr] for tr in wrap.iter('tr')]
    if tuple(rows[0]) != COLUMNS:
        raise ValueError(f'Unexpected header of the integrin table: {rows[0]!r}')
    seen: set[str] = set()
    for cells in rows[1:]:
        if len(cells) != len(COLUMNS):
            raise ValueError(f'Row of the integrin table with {len(cells)} cells: {cells!r}')
        row = dict(zip(COLUMNS, cells))
        if row['Gene symbol'] in seen:
            raise ValueError(f'Subunit {row["Gene symbol"]!r} occurs twice')
        seen.add(row['Gene symbol'])
        if not row['Gene symbol'] or not row['Protein accession number']:
            raise ValueError(f'Subunit without symbol or accession: {cells!r}')
        yield {
            'symbol': row['Gene symbol'],
            'name': row['Protein name'],
            'synonyms': row['Synonyms'] or None,
            'gene_accession': row['Gene accession number'],
            'location': row['Location (chromosome)'],
            'accession': row['Protein accession number'],
        }
    if not seen:
        raise ValueError('The integrin table has no rows')


def iter_family(opener, **kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """The family: one record with the ``members`` (the subunits, see
    :func:`iter_subunits`)."""
    yield {'members': list(iter_subunits(opener, **kwargs))}
