"""
Reader for the HIPPIE interaction table (``hippie_current.txt``).

HIPPIE (Human Integrated Protein-Protein Interaction rEference) distributes
its interactions as a tab separated table with six columns and no header
line, one row per pair of UniProt entries:

    AL1A1_HUMAN  216  AL1A1_HUMAN  216  0.76  experiments:in vivo,Two-hybrid;pmids:12081471,16189514;sources:HPRD,BioGRID

The columns are the UniProt entry name and the Entrez Gene ID of each
protein, the confidence score and the evidence of the pair, packed into one
cell as ``key:value,value;key:value,value``. The keys are ``experiments``
(names of detection techniques, free text of the source), ``pmids``,
``sources`` (the databases and studies the pair comes from) and ``species``
(other organisms in which the interaction was found; only some rows have it).
An entry name cell can hold several UniProt entries of one gene, separated by
a comma, and can be empty.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Generator
import logging
import re
from typing import Any

from pypath.inputs_v2.parsers.base import iter_tsv

__all__ = ['COLUMNS', 'DETAIL_KEYS', 'iter_hippie', 'parse_details']

_log = logging.getLogger(__name__)

# Columns of the table, in file order. The file has no header line.
COLUMNS = (
    'entry_names_a',
    'entrez_a',
    'entry_names_b',
    'entrez_b',
    'score',
    'details',
)
# The keys of the packed last column, as the source writes them.
DETAIL_KEYS = ('experiments', 'pmids', 'sources', 'species')

_ENTRY_NAME_RE = re.compile(r'[A-Z0-9]+_[A-Z0-9]+')
_SCORE_RE = re.compile(r'(?:0(?:\.\d+)?|1(?:\.0+)?)')


def _entry_names(cell: str) -> list[str]:
    """The UniProt entry names of a cell, in order and without repeats."""
    names = []
    for name in (part.strip() for part in cell.split(',')):
        if not name:
            continue
        if not _ENTRY_NAME_RE.fullmatch(name):
            raise ValueError(f'HIPPIE: not a UniProt entry name: {name!r}')
        if name not in names:
            names.append(name)
    return names


def parse_details(cell: str) -> dict[str, list[str]]:
    """Unpack the last column into the lists of ``DETAIL_KEYS``.

    The lists are as written, in file order; a key the cell lacks gives an
    empty list. An empty cell gives four empty lists. Any other key, or a key
    given twice, is an error: the layout of the file changed.
    """
    details: dict[str, list[str]] = {key: [] for key in DETAIL_KEYS}
    seen = set()
    for part in cell.split(';') if cell else []:
        key, colon, values = part.partition(':')
        if not colon or key not in details or key in seen:
            raise ValueError(f'HIPPIE: unexpected evidence item {part!r}')
        seen.add(key)
        details[key] = [value.strip() for value in values.split(',') if value.strip()]
    return details


def iter_hippie(
    opener,
    stats: Counter | None = None,
    **_kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """Yield one dictionary per row of the table.

    The fields are:

    - ``entry_names_a``, ``entry_names_b``: lists of UniProt entry names (can
      be empty);
    - ``entrez_a``, ``entrez_b``: Entrez Gene IDs, as text;
    - ``score``: the confidence score as written, a number from 0 to 1;
    - ``experiments``, ``pmids``, ``sources``, ``species``: lists, as written
      (see :func:`parse_details`).

    A row with another number of fields, an Entrez ID that is not a number, a
    score outside 0 to 1, a name that is not a UniProt entry name, or an
    unknown evidence key is an error, and so is a table with no rows. ``stats``,
    if given, is filled with the counts ``rows`` and ``without_evidence``
    (rows whose last column is empty), which are also logged.
    """
    stats = Counter() if stats is None else stats
    for row in iter_tsv(opener, fieldnames=COLUMNS, strict=True, quotechar=None):
        stats['rows'] += 1
        for column in ('entrez_a', 'entrez_b'):
            if not row[column].isdigit():
                raise ValueError(f'HIPPIE: not an Entrez Gene ID: {row[column]!r}')
        if not _SCORE_RE.fullmatch(row['score']):
            raise ValueError(f'HIPPIE: score not in 0 to 1: {row["score"]!r}')
        details = parse_details(row['details'])
        if not row['details']:
            stats['without_evidence'] += 1
        yield {
            'entry_names_a': _entry_names(row['entry_names_a']),
            'entrez_a': row['entrez_a'],
            'entry_names_b': _entry_names(row['entry_names_b']),
            'entrez_b': row['entrez_b'],
            'score': row['score'],
            **details,
        }
    if not stats['rows']:
        raise ValueError('HIPPIE table has no rows')
    _log.info('HIPPIE table: %s', dict(sorted(stats.items())))
