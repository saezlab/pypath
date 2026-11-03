"""
Reader for the DoRothEA ``entire_database`` table.

The table is an R data frame in the ``entire_database.rda`` file of the
DoRothEA R package, read by the shared R data reader. It holds every human
transcription factor - target gene pair of DoRothEA, with a confidence level
from A (best supported) to E; :func:`iter_level` selects the rows of one
level.
"""

from __future__ import annotations

from collections.abc import Generator
from typing import Any

from pypath.inputs_v2.parsers.base import iter_rdata

__all__ = ['CONFIDENCE_LEVELS', 'iter_level']

# Confidence levels of the table, from the best supported to computational
# predictions only.
CONFIDENCE_LEVELS = ('A', 'B', 'C', 'D', 'E')


def iter_level(
    opener,
    level: str,
    **kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """Rows of the ``entire_database`` data frame with one confidence level.

    The fields of a row are the columns of the data frame: ``tf`` and
    ``target`` (HGNC symbols), ``mor`` (mode of regulation: 1, -1 or 0 for
    not known), ``confidence`` (the level), ``is_evidence_curated``,
    ``is_evidence_chip_seq``, ``is_evidence_tfbs`` and ``is_evidence_inferred``
    (booleans), ``which_curated``, ``which_chip_seq``, ``which_tfbs`` and
    ``which_inferred`` (the data sources of each line of evidence, joined by
    commas; ``none`` for no source) and ``pubmed_id`` (references).

    A row with a level outside A to E belongs to no dataset, so it is an error
    here rather than a row that quietly disappears.
    """
    for row in iter_rdata(opener, name='entire_database', **kwargs):
        row_level = row.get('confidence')
        if row_level not in CONFIDENCE_LEVELS:
            raise ValueError(f'Unknown DoRothEA confidence level: {row_level!r}')
        if row_level == level:
            yield row
