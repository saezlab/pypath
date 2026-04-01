"""
Common raw data parsers for simple file formats.

These parsers handle standard formats like CSV, TSV, JSON, JSONL, and SQLite.
"""

from __future__ import annotations

from collections.abc import Generator
import csv
import itertools
import json
from pathlib import Path
import sqlite3
from typing import Any
import os


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


def iter_sqlite(
    opener,
    table_name: str | None = None,
    sqlite_path: Path | None = None,
    db_rel_path: str | None = None,
    query: str | None = None,
    **_kwargs: Any
) -> Generator[dict[str, Any], None, None]:
    """
    Iterates over a table or query in a SQLite database.
    If the database is within an archive (e.g., .tar.gz), it extracts it
    to a local cache path first.

    Args:
        opener: The opener object for handling database connections.
        table_name: The name of the table to retrieve data from (if query not provided).
        sqlite_path: The local path to the cached SQLite database.
        db_rel_path: Optional relative path to the database file within the archive.
        query: Optional custom SQL query to execute.
        max_records: Optional maximum number of records to retrieve.

    Yields:
        Dictionary representing a row from the table or query.
    """

    if sqlite_path and not os.path.exists(sqlite_path):

        if not opener or not opener.result:

            return

        _extract_sqlite_from_opener(opener, sqlite_path, db_rel_path)

    results = _iter_table(sqlite_path, table_name, query=query)

    yield from results


def _extract_sqlite_from_opener(opener, target_path: Path, db_rel_path: str | None) -> None:
    """Extracts a SQLite database file from the opener to a local path."""
    if isinstance(opener.result, dict):
        if not db_rel_path:
            raise ValueError("Relative path for file in archive is required for extracted archive results.")

        if db_rel_path not in opener.result:
            raise FileNotFoundError(f"Could not find {db_rel_path} in archive.")

        db_handle = opener.result[db_rel_path]
    else:
        db_handle = opener.result

    print(f"Extracting SQLite database to: {target_path}")
    target_path.parent.mkdir(parents=True, exist_ok=True)

    with open(target_path, 'wb') as fp:
        while chunk := db_handle.read(1024 ** 2):
            fp.write(chunk)


def _iter_table(
    sqlite_path: Path,
    table_name: str | None = None,
    query: str | None = None,
) -> Generator[dict[str, Any], None, None]:
    """Pure generator for yielding rows from a SQLite table or query."""
    connection = None
    try:
        connection = sqlite3.connect(sqlite_path)
        connection.row_factory = sqlite3.Row
        cursor = connection.cursor()

        if not query:
            if not table_name:
                raise ValueError("Either table_name or query must be provided.")
            query = f"SELECT * FROM {table_name};"

        cursor.execute(query)

        for row in cursor:
            yield dict(row)

    except sqlite3.Error as e:
        print(f"SQLite error during iteration: {e}")
        raise
    finally:
        if connection:
            connection.close()
