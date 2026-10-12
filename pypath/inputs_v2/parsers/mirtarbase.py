"""
Read the miRTarBase download tables.

Source files (comma separated, UTF-8 with a byte order mark, one experimental
support block of a miRNA - target gene interaction per line):

    miRTarBase_SE_WR.csv    strong experimental evidence (reporter assay or Western blot)
    miRTarBase_WE_Other.csv weak experimental evidence other than CLIP-Seq

    miRTarBase ID,miRNA,Species (miRNA),Target Gene,Target Gene (Entrez ID),Species (Target Gene),Experiments,Support Type,References (PMID)
    MIRT003105,hsa-miR-122-5p,hsa,SLC7A1,6541.0,hsa,Luciferase reporter assay//Western blot,Functional MTI,17179747.0

The Entrez Gene ID and the PubMed ID are written as floats (``6541.0``). The
``Experiments`` cell joins the experiment labels of the block with ``//``.
"""

from __future__ import annotations

from collections.abc import Generator
import re
from typing import Any

from pypath.inputs_v2.parsers.base import iter_csv

# Header of the file, in source order.
COLUMNS = (
    'miRTarBase ID',
    'miRNA',
    'Species (miRNA)',
    'Target Gene',
    'Target Gene (Entrez ID)',
    'Species (Target Gene)',
    'Experiments',
    'Support Type',
    'References (PMID)',
)

_FLOAT_ID = re.compile(r'^(\d+)(?:\.0+)?$')


def _clean(value: Any) -> str:
    return (value or '').strip()


def _integer(value: Any, column: str, number: int) -> str | None:
    """An ID the file writes as a float (``6541.0``), as digits; None if empty."""
    text = _clean(value)
    if not text:
        return None
    match = _FLOAT_ID.match(text)
    if not match:
        raise ValueError(f'Row {number}: {column} is not a number: {text!r}')
    return match.group(1)


def iter_interactions(
    opener: Any,
    **kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """
    Yield the support blocks of a miRTarBase table.

    Keys: ``mirtarbase_id``, ``mirna`` (the mature miRNA name as written, with
    surrounding spaces removed), ``mirna_species`` and ``target_species`` (the
    file's three or four letter organism codes), ``target_gene`` (symbol),
    ``entrez_id`` (digits, None when empty), ``experiments`` (the labels of the
    ``Experiments`` cell as a list, empty labels left out), ``support_type`` and
    ``pmid`` (digits as written; the file has ``0`` and ``1`` for some blocks
    without a publication, kept as they are). Every line is one block; the same
    block is repeated in the file where the source lists it more than once.

    The header must have exactly the columns of the source; an empty table is an
    error.
    """
    count = 0
    for number, row in enumerate(iter_csv(opener, **kwargs), start=2):
        row = {
            (name or '').lstrip('﻿').strip(): value
            for name, value in row.items()
        }
        if set(row) != set(COLUMNS):
            raise ValueError(
                f'Unexpected miRTarBase columns: {sorted(map(str, row))}'
            )
        if any(row[name] is None for name in COLUMNS):
            raise ValueError(f'Row {number} has too few fields')
        count += 1
        yield {
            'mirtarbase_id': _clean(row['miRTarBase ID']),
            'mirna': _clean(row['miRNA']),
            'mirna_species': _clean(row['Species (miRNA)']),
            'target_gene': _clean(row['Target Gene']),
            'entrez_id': _integer(row['Target Gene (Entrez ID)'], 'Entrez ID', number),
            'target_species': _clean(row['Species (Target Gene)']),
            'experiments': [
                label.strip()
                for label in row['Experiments'].split('//')
                if label.strip()
            ],
            'support_type': _clean(row['Support Type']),
            'pmid': _integer(row['References (PMID)'], 'PMID', number),
        }
    if not count:
        raise ValueError('miRTarBase table has no rows')
