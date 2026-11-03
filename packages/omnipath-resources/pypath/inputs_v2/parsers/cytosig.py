"""
Parser of the CytoSig signature matrix.

``signature.centroid`` is a tab separated matrix of 4,881 genes (rows) by 43
cytokines (columns). Its header row names only the 43 cytokines and not the
first column, which holds the gene symbols (the file is read with ``index_col=0``
by CytoSig itself), so a row has one field more than the header.
"""

from __future__ import annotations

from collections.abc import Generator
import io
from types import SimpleNamespace
from typing import Any

from pypath.inputs_v2.parsers.base import iter_tsv, read_text

GENE = 'gene'


def iter_signature(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """One row for each cell of the matrix: ``gene``, ``cytokine`` (the column
    name as written) and ``value`` (the cell, as written).

    A row whose number of fields is not that of the header plus one, a header
    that repeats a name, an empty cell and a table without rows raise.
    """
    text = read_text(opener)
    header, _, _ = text.partition('\n')
    cytokines = header.rstrip('\r').split('\t')
    if not header.strip() or len(cytokines) < 2:
        raise ValueError('CytoSig signature matrix has no header')
    rows = iter_tsv(
        SimpleNamespace(result=io.StringIO(text)),
        fieldnames=(GENE, *cytokines),
        skip_rows=1,
        strict=True,
        quotechar=None,
        required=(GENE, *cytokines),
    )
    for row in rows:
        for cytokine in cytokines:
            value = row[cytokine]
            if value is None or not value.strip():
                raise ValueError(f'CytoSig signature has an empty cell: {row[GENE]} / {cytokine}')
            yield {'gene': row[GENE], 'cytokine': cytokine, 'value': value}


def iter_cytokines(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """The cytokines (column names) of the matrix, as written, in column order."""
    header, _, _ = read_text(opener).partition('\n')
    cytokines = header.rstrip('\r').split('\t')
    if not header.strip() or len(cytokines) < 2 or len(set(cytokines)) != len(cytokines):
        raise ValueError('CytoSig signature matrix has no usable header')
    for cytokine in cytokines:
        yield {'cytokine': cytokine}
