"""
Common raw data parsers for simple file formats.

These parsers handle standard formats like CSV, TSV, JSON, and JSONL.
"""

from __future__ import annotations

from collections.abc import Generator
import csv
import json
from typing import Any


def _first_handle(opener) -> Any | None:
    """Extract the first file handle from an opener result."""
    if not opener or not opener.result:
        return None
    if isinstance(opener.result, dict):
        return next(iter(opener.result.values()), None)
    return opener.result


def _raw(opener, delimiter: str = ',', **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """
    Parse CSV files with configurable delimiter.

    Args:
        opener: File opener from download_and_open
        delimiter: Field delimiter (default: ',')

    Yields:
        Dictionary for each row
    """
    handle = _first_handle(opener)
    if not handle:
        return
    yield from csv.DictReader(handle, delimiter=delimiter)


def iter_csv(opener, delimiter: str = ',', **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Parse CSV files with configurable delimiter."""
    yield from _raw(opener, delimiter=delimiter, **_kwargs)


def iter_tsv(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Parse TSV (tab-separated values) files."""
    yield from _raw(opener, delimiter='\t', **_kwargs)


def iter_semicolon(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Parse semicolon-delimited files."""
    yield from _raw(opener, delimiter=';', **_kwargs)


def iter_json(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """
    Parse JSON files.

    Yields individual records if the JSON is a list, otherwise yields the whole object.
    """
    handle = _first_handle(opener)
    if not handle:
        return
    data = json.load(handle)
    if isinstance(data, list):
        yield from data
    else:
        yield data


def iter_jsonl(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Parse JSONL (JSON Lines) files - one JSON object per line."""
    handle = _first_handle(opener)
    if not handle:
        return
    for line in handle:
        if line.strip():
            yield json.loads(line)
