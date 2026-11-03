"""
Reader for the STRING protein links files.

The ``links.detailed`` files of STRING (functional, physical and regulatory
network) are space separated tables with a header, one protein pair per row:
the two STRING identifiers, the scores of the evidence channels, and the
combined score. Scores are probabilities multiplied by 1000.
"""

from __future__ import annotations

from collections.abc import Generator
from typing import Any

from pypath.inputs_v2.parsers.base import iter_csv

__all__ = ['iter_links']


def _is_score(value: object) -> bool:
    """A whole number from 0 to 1000 written in ASCII digits."""
    return (
        isinstance(value, str)
        and value.isascii()
        and value.isdigit()
        and int(value) <= 1000
    )


def iter_links(
    opener,
    columns: tuple[str, ...] = (),
    minimum_score: int = 0,
    directed: bool = False,
    **_kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """Yield the rows of a STRING links file with at least ``minimum_score``.

    Fields of a row, all strings as written: ``protein1`` and ``protein2``
    (STRING identifiers, the NCBI taxon and the protein identifier of its
    source genome joined by a dot: ``9606.ENSP00000000233``), one field for
    each name in ``columns`` (an evidence channel score, or for a regulatory
    file the score of positive or negative regulation) and
    ``combined_score``. The file must have these columns.

    The files of the functional and the physical network list every pair
    twice, as A B and as B A, with the same scores. One row per pair is
    yielded, the one whose first identifier sorts first, and a file in which
    the two orientations do not match in number is rejected. A regulatory
    file is ``directed``: a row is the link from the first protein to the
    second, and every row is kept.

    Every score of every row (the channels in ``columns`` and the combined
    score) must be a whole number from 0 to 1000; any other value raises
    ``ValueError``, so that no score is dropped later for being malformed.
    """
    scores = (*columns, 'combined_score')
    required = ('protein1', 'protein2', *scores)
    rows = forward = mirrored = 0
    for row in iter_csv(opener, delimiter=' '):
        if not rows:
            missing = [column for column in required if column not in row]
            if missing:
                raise ValueError(
                    f'STRING links file lacks expected columns: {", ".join(missing)}'
                )
        rows += 1
        for column in scores:
            if not _is_score(row[column]):
                raise ValueError(
                    f'STRING links file row {rows} has {row[column]!r} as the '
                    f'{column} score, not a whole number from 0 to 1000'
                )
        if not directed and row['protein1'] != row['protein2']:
            if row['protein1'] > row['protein2']:
                mirrored += 1
                continue
            forward += 1
        if int(row['combined_score']) >= minimum_score:
            yield row
    if not rows:
        raise ValueError('STRING links file has no rows')
    if forward != mirrored:
        raise ValueError(
            'STRING links file is not symmetric: '
            f'{forward} pairs in one orientation, {mirrored} in the other'
        )
