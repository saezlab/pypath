"""
Parser for MRClinksDB metabolite-receptor interaction data.

Reads a tab-separated MRClinksDB file and yields rows as dicts with
normalised snake_case keys, ready for FieldConfig transformations.
"""

from __future__ import annotations

import csv
import re
from collections.abc import Generator
from typing import Any

from pypath.inputs_v2.base import _first_handle


def _normalise_header(raw: str) -> str:
    """Convert a raw header string to a lowercase snake_case key."""
    return re.sub(r'[^a-z0-9]+', '_', raw.strip().lower()).strip('_')


def iter_mrclinksdb_interactions(
        opener: Any,
        **_kwargs: Any,
) -> Generator[dict[str, str], None, None]:
    """
    Parse a MRClinksDB interaction TSV and yield normalised row dicts.

    Like ``iter_tsv``, this function uses the file's own header row as
    dict keys, so column-order differences between organism files (e.g.
    the mouse file swaps the gene-ID and UniProt-ID columns) are handled
    transparently — no positional swap detection is required.

    The difference from ``iter_tsv`` is that header strings are normalised
    to snake_case before being used as keys (e.g. ``"Receptor UniProt ID"``
    → ``"receptor_uniprot_id"``, ``"PubChem CID/SID"`` → ``"pubchem_cid_sid"``).
    This keeps FieldConfig selectors in the schema module stable and
    readable regardless of how the source file capitalises its headers.

    Multi-value splitting (PubChem CID extraction, UniProt complex
    splitting, PMID lists) is left to FieldConfig transforms in the
    schema module.

    Args:
        opener: File opener returned by ``Download.open()``.
        **_kwargs: Ignored; accepted for interface compatibility.

    Yields:
        Dict mapping normalised column name → raw cell string.
    """
    handle = _first_handle(opener)
    if not handle:
        return

    reader = csv.reader(handle, delimiter='\t')
    raw_header = next(reader, None)
    if raw_header is None:
        return

    fieldnames = [_normalise_header(h) for h in raw_header]

    for fields in reader:
        if not any(fields):
            continue
        yield dict(zip(fieldnames, fields))
