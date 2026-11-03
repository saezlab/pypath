"""
Read the miR2Disease table of experimentally verified miRNA - target pairs.

Source file (``miRtar.txt``, tab separated, ISO-8859-1, CRLF line ends):

    The following miRNA-target pairs, validated by luciferase reporter experiments, were extracted manually from the corresponding reference. A few target genes not documented here may be available in "TarBase" and "miRecords" database system.

    miRNA   Validated target        Pub Date        Reference
    hsa-let-7a      KRAS    2008    A SNP in a let-7 microRNA complementary site in the KRAS 3' untranslated region increases non-small cell lung cancer risk.

Two lines of text and a blank line come before the header. The reference is the
title of the paper; the file gives no PubMed ID. A title that was broken over
several lines in the source continues on lines without tabs, and a few titles are
wrapped in double quotes.
"""

from __future__ import annotations

from collections.abc import Generator
import re
from typing import Any

from pypath.inputs_v2.parsers.base import iter_tsv

COLUMNS = ('miRNA', 'Validated target', 'Pub Date', 'Reference')

_YEAR = re.compile(r'^\d{4}$')


def _record(row: dict[str, Any]) -> dict[str, Any]:
    mirna = (row['miRNA'] or '').strip()
    species, dash, _ = mirna.partition('-')
    if not mirna or not dash:
        raise ValueError(f'miR2Disease miRNA name without organism code: {mirna!r}')
    year = (row['Pub Date'] or '').strip()
    if not _YEAR.match(year):
        raise ValueError(f'miR2Disease publication year is not a year: {year!r}')
    reference = (row['Reference'] or '').strip()
    if len(reference) > 1 and reference.startswith('"') and reference.endswith('"'):
        reference = reference[1:-1].strip()
    return {
        'mirna': mirna,
        'mirna_species': species.lower(),
        'target': (row['Validated target'] or '').strip(),
        'year': year,
        'reference': reference,
    }


def iter_targets(
    opener: Any,
    **kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """
    Yield the pairs of the table.

    Keys: ``mirna`` (the miRNA name, spaces around it removed), ``mirna_species``
    (the organism code before the first ``-`` of the name, lowercase),
    ``target`` (the name of the target as written, spaces around it removed; it
    is a gene or protein name, not always a gene symbol), ``year`` (four digits)
    and ``reference`` (the paper title on one line, without the double quotes
    that wrap some titles).

    The three lines before the header are skipped and the header must have the
    four columns of the source. A line that does not start a pair continues the
    title of the pair before it. An empty table is an error.
    """
    pending = None
    count = 0
    # No quote character: a quote opens a title in the source but is not always
    # closed on the same line, and titles contain quotes of their own.
    for row in iter_tsv(opener, skip_rows=2, quotechar=None, **kwargs):
        if tuple(row) != COLUMNS:
            raise ValueError(f'Unexpected miR2Disease columns: {list(row)}')
        if row['Validated target'] is None:
            # Continuation of a title that the source broke over lines.
            if pending is None:
                raise ValueError('miR2Disease table starts with a continuation line')
            pending['Reference'] = f"{pending['Reference']} {row['miRNA']}"
            continue
        if pending is not None:
            yield _record(pending)
            count += 1
        pending = row
    if pending is not None:
        yield _record(pending)
        count += 1
    if not count:
        raise ValueError('miR2Disease table has no rows')
