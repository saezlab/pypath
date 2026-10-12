"""
Parser of the literature curated TransmiR file.

``db_literaturev3.txt`` of TransmiR v3.0 is tab separated, without a header,
one transcription factor - microRNA regulation per line, with nine columns:

    AGL15  csi-MIR156a  (empty)  (empty)  Activation  35661206  literature  embryogenic callus  C.sinensis

 1. ``tf``: the name of the transcription factor, as the curator wrote it
 2. ``mirna``: the name of the microRNA, mostly a miRBase precursor name
    (``hsa-mir-21``), written ``MIR156a`` for plants
 3. ``tss``: transcription start site of the microRNA gene, empty in this file
 4. ``site``: transcription factor binding site, empty in this file
 5. ``action``: the action type (``Activation``, ``Repression``,
    ``Regulation``), with a remark in brackets for some lines
 6. ``pmid``: the PubMed ID of the publication
 7. ``evidence``: ``literature``
 8. ``tissue``: the tissue or cell line, empty for most lines
 9. ``species``: the species of the regulation, as TransmiR writes it

A few cells have a space at the end; it is removed. Lines are not
deduplicated.
"""

from __future__ import annotations

from collections.abc import Generator
from typing import Any

from pypath.inputs_v2.parsers.base import iter_tsv

COLUMNS = (
    'tf',
    'mirna',
    'tss',
    'site',
    'action',
    'pmid',
    'evidence',
    'tissue',
    'species',
)


def iter_literature(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """
    One dictionary per line of the literature file, keyed by ``COLUMNS``.

    A line that does not have nine columns, or a file with no line, is an
    error.
    """
    count = 0
    for row in iter_tsv(
        opener,
        fieldnames=COLUMNS,
        strict=True,
        quotechar=None,
    ):
        count += 1
        yield {key: (value or '').strip() for key, value in row.items()}
    if not count:
        raise ValueError('TransmiR literature file has no lines')
