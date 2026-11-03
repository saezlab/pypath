"""
Read the TRRUST regulatory relationship files.

Source files (tab separated, no header; one file per species):

    trrust_rawdata.human.tsv
    trrust_rawdata.mouse.tsv

    AATF    BAX     Repression  22909821
    ABL1    BAX     Activation  11753601
    AATF    CDKN1A  Unknown     17157788

The columns are the transcription factor, the target gene (both official gene
symbols of the file's species), the mode of regulation and the PubMed IDs of
the articles the relationship was curated from, joined by ``;``.
"""

from __future__ import annotations

from collections.abc import Generator
from typing import Any

from pypath.inputs_v2.parsers.base import iter_tsv

# The files have no header row; these names are ours, in source order.
COLUMNS = ('tf', 'target', 'mode', 'pmids')


def iter_relationships(
    opener,
    *,
    taxon_id: str,
    **kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """
    Yield the lines of a TRRUST file with the file's species added.

    Keys: ``tf``, ``target``, ``mode`` and ``pmids`` as the file has them, and
    ``taxon_id``, the NCBI Taxonomy ID of the file. A line with another number
    of fields than four is an error, and so is a file without lines.
    """
    count = 0
    for row in iter_tsv(
        opener,
        fieldnames=COLUMNS,
        strict=True,
        quotechar=None,
        **kwargs,
    ):
        row['taxon_id'] = taxon_id
        count += 1
        yield row
    if not count:
        raise ValueError('TRRUST file has no rows')
