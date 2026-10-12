"""
Read the GFF file of the ABS database.

Source file (``abs.gff``, tab separated, no header; one annotated binding site
of a transcription factor in a promoter sequence per line, with a tab ending
every line):

    M21390     A0012   P53   326   335   .   +   .   # MM:  ATACAAGGCC
    extracted04  A0001   SRF   262   271   .   +   .   # HS:  TCCTTCTTTG

The columns are the promoter sequence (a GenBank accession, or ``extractedNN``
for a promoter that the authors extracted from a genome sequence), the ABS
entry (``A0012``, a set of orthologous promoters), the name of the
transcription factor, the first and last position of the site in the 500
nucleotides of the promoter, the score and the frame (always ``.``), the
strand, and a comment with the species of the promoter (``HS``, ``MM``,
``RN`` or ``GG``) and the sequence of the site.
"""

from __future__ import annotations

from collections.abc import Generator
import re
from typing import Any

from pypath.inputs_v2.parsers.base import iter_tsv

COLUMNS = (
    'sequence',
    'entry',
    'factor',
    'start',
    'end',
    'score',
    'strand',
    'frame',
    'comment',
    'trailing',
)
_COMMENT = re.compile(r'# ([A-Z]{2}):\s+([ACGTN]+)')


def iter_sites(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """
    Yield the binding sites of the GFF file.

    Keys: ``sequence``, ``entry`` and ``factor`` (as written), ``start`` and
    ``end`` (integers, one-based, in the promoter sequence), ``strand`` (``+``
    or ``-``), ``species`` (the code of the comment) and ``site`` (the
    sequence of the site). A line with another number of fields, a score or
    frame other than ``.``, a comment that is not a species and a sequence, a
    position that is not a number or whose end precedes its start, and a file
    without lines are errors.
    """
    for number, row in enumerate(
        iter_tsv(opener, fieldnames=COLUMNS, strict=True, quotechar=None, required=COLUMNS),
        1,
    ):
        comment = _COMMENT.fullmatch(row['comment'].strip())
        if not comment or row['score'] != '.' or row['frame'] != '.':
            raise ValueError(f'ABS line {number} is not a binding site: {row}')
        if row['strand'] not in ('+', '-'):
            raise ValueError(f'ABS line {number} has the strand {row["strand"]!r}')
        start, end = int(row['start']), int(row['end'])
        if not 0 < start <= end:
            raise ValueError(f'ABS line {number} has the positions {start}-{end}')
        yield {
            'sequence': row['sequence'],
            'entry': row['entry'],
            'factor': row['factor'],
            'start': start,
            'end': end,
            'strand': row['strand'],
            'species': comment[1],
            'site': comment[2],
        }
