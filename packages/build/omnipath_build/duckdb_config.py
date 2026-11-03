"""Memory and thread budgets for resource-build DuckDB connections."""

import os
import re


def build_memory_limit(explicit: str | None = None) -> str:
    value = (
        explicit if explicit is not None else os.environ.get("OMNIPATH_BUILD_DUCKDB_MEMORY", "1GB")
    )
    match = re.fullmatch(r"\s*(\d+(?:\.\d*)?|\.\d+)\s*(B|[KMGTPE]I?B)\s*", value, re.IGNORECASE)
    if not match or float(match[1]) <= 0:
        raise ValueError("Build DuckDB memory must be a positive size, for example 1GB or 512MiB")
    return match[1] + match[2].upper()


def configure_memory(connection, limit: str | None = None) -> None:
    connection.execute("SET memory_limit = ?", [build_memory_limit(limit)])


GIB = 1024**3


def finalize_budget(ram_bytes: int, cpus: int) -> tuple[str, int]:
    """DuckDB memory limit and threads for the finalizer inside a hard memory reservation.

    DuckDB's memory_limit covers its buffer pool only. Per-thread operator state, the Python
    process and the page cache come on top, so a limit near the reservation gets the job
    killed instead of spilling. Half the reservation for DuckDB and about one thread per GiB
    keep the total under the kill threshold; DuckDB spills the rest to its temp directory.
    """
    limit = max(256 * 1024**2, ram_bytes // 2)
    threads = max(1, min(cpus, ram_bytes // GIB))
    return f"{limit}B", threads


def build_threads() -> int:
    value = int(os.environ.get("OMNIPATH_BUILD_DUCKDB_THREADS", "4"))
    if value < 1:
        raise ValueError("Build DuckDB threads must be positive")
    return value
