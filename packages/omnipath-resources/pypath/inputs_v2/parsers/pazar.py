"""
Parser of the PAZAR transcription factor - target gene lists.

The archive ``tftargets.zip`` of the PAZAR site holds one tab separated file
per project, without a header (the project ``ABS`` is ``pazar_ABS_20120522.csv``;
an empty project is an empty file). Every line is one assay of one transcription
factor on one gene sequence, with twelve columns:

    TF0000324  ENST00000265354  SRF_HUMAN  GS0000616  ENSG00000143632  1  229566992  229569845  Homo sapiens  ABS  12760745  NA

 1. ``tf_id``: the PAZAR transcription factor record (``TF`` and seven digits)
 2. ``tf_transcript``: the Ensembl transcript of the factor, or in a few lines an
    Ensembl gene, a FlyBase gene or transcript, or a worm sequence name
 3. ``tf_name``: the name of the factor, mostly UniProt style (``SRF_HUMAN``)
 4. ``gene_sequence_id``: the PAZAR gene sequence record (``GS`` and seven digits)
 5. ``gene``: the Ensembl gene of the target, or a FlyBase gene or a worm
    sequence name
 6. ``chromosome``: of the gene sequence, ``-`` when it has no coordinates
 7. ``start`` and 8. ``end``: coordinates of the gene sequence
 9. ``species``: the species of the target gene
10. ``project``: the PAZAR project of the line
11. ``pmid``: the PubMed ID of the publication, empty in a few lines
12. ``method``: the assay, as ``class::method``; ``NA`` where the project gives
    none

Two repairs are made: surrounding spaces are removed, and ``a :: b`` is written
``a::b`` (one project writes the separator with spaces). FlyBase identifiers
that are in upper case (``FBGN0039044``) are written as FlyBase writes them
(``FBgn0039044``).

Lines are not deduplicated: a line repeated in a project is a repeated
observation.
"""

from __future__ import annotations

from collections.abc import Generator
import re
from types import SimpleNamespace
from typing import Any

from pypath.inputs_v2.parsers.base import iter_tsv

COLUMNS = (
    'tf_id',
    'tf_transcript',
    'tf_name',
    'gene_sequence_id',
    'gene',
    'chromosome',
    'start',
    'end',
    'species',
    'project',
    'pmid',
    'method',
)
_FLYBASE = re.compile(r'^FB(GN|TR)(\d+)$')


def _flybase(value: str) -> str:
    match = _FLYBASE.fullmatch(value)
    if match:
        return f"FB{match.group(1).lower()}{match.group(2)}"
    return value


def iter_tf_targets(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """
    One dictionary per line of every file of the archive.

    The keys are ``COLUMNS``. Files are read in the order of the archive. An
    archive without a line, or a member that is not a table of twelve columns,
    is an error.
    """
    if not opener or not opener.result or not isinstance(opener.result, dict):
        raise ValueError('PAZAR download is not an archive of files')
    count = 0
    for name, handle in opener.result.items():
        for row in iter_tsv(
            SimpleNamespace(result=handle),
            fieldnames=COLUMNS,
            strict=True,
            quotechar=None,
        ):
            record = {key: (value or '').strip() for key, value in row.items()}
            record['method'] = re.sub(r'\s*::\s*', '::', record['method'])
            record['tf_transcript'] = _flybase(record['tf_transcript'])
            record['gene'] = _flybase(record['gene'])
            count += 1
            yield record
    if not count:
        raise ValueError('PAZAR archive has no lines')
