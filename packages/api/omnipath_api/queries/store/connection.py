"""DuckDB connection management and read_parquet expression builders."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Sequence
import duckdb


def get_connection(memory_limit: str = "4GB", *, threads: int = 2) -> duckdb.DuckDBPyConnection:
    """Create configured in-memory DuckDB connection."""
    threads = int(os.getenv("OMNIPATH_DUCKDB_THREADS", str(threads)))
    if threads < 1:
        raise ValueError("OMNIPATH_DUCKDB_THREADS must be positive")
    con = duckdb.connect(database=":memory:")
    try:
        con.execute(
            "SET memory_limit = ?", [os.getenv("OMNIPATH_DUCKDB_MEMORY_LIMIT", memory_limit)]
        )
        con.execute("SET threads = ?", [threads])
        con.execute("SET arrow_large_buffer_size=true")
        # Footers of the immutable resource files are read once, not on every query.
        con.execute("SET parquet_metadata_cache=true")
    except Exception:
        con.close()
        raise
    return con


def format_read_parquet(
    paths: Sequence[str | Path],
    union_by_name: bool = True,
    *,
    filename: bool = False,
    file_row_number: bool = False,
) -> str:
    if not paths:
        return "read_parquet([], union_by_name=true)"
    quoted = ", ".join(sql_literal(str(path)) for path in paths)
    union_flag = ", union_by_name=true" if union_by_name else ""
    options = (
        union_flag
        + (", filename=true" if filename else "")
        + (", file_row_number=true" if file_row_number else "")
    )
    return f"read_parquet([{quoted}]{options})"


def sql_literal(value: str) -> str:
    """Escape the rare SQL literal that cannot use a bound query parameter."""
    return "'" + value.replace("'", "''") + "'"
