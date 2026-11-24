"""The query budget of an API process: how many queries run at once, and the DuckDB
threads and memory of the database they share, from the CPUs and memory its container
gets.

One query slot per CPU, with one thread each: the explorer's queries are small, so
waiting for a slot costs more than a second thread per query gains. In load tests on
7 CPUs and 13 GB (the explorer's extended journey, 64 people), six slots of one
thread answered the slowest 5% of steps in 0.31 s, three slots of two in 1.07 s, at the
same CPU per request; 13 GB and 32 GB answered alike. DuckDB gets most of the
container's memory, the rest is for Python, Arrow and response buffers.

Each part can be set: OMNIPATH_QUERY_CONCURRENCY (slots), OMNIPATH_DUCKDB_THREADS
(threads per slot) and OMNIPATH_DUCKDB_MEMORY_LIMIT (memory per slot). An empty setting
counts as unset."""

from __future__ import annotations

from dataclasses import dataclass
import os
from pathlib import Path
import re

from omnipath_client.budget import cgroup_memory_limit, cpu_count

# The share of the container's memory DuckDB may use (tested: the whole process peaked at
# 9 GB of 13 GB with DuckDB allowed 12 GB).
MEMORY_SHARE = 0.75
# DuckDB's memory per slot without a container memory limit.
FALLBACK_MEMORY = "1GB"

_SIZE = re.compile(r"\s*(\d+(?:\.\d+)?)\s*([KMGT]?)(i?)B?\s*", re.I)


def scaled(size: str, factor: int) -> str:
    """A DuckDB memory size times ``factor`` (``3500MB`` times 3 is ``10500MB``)."""
    match = _SIZE.fullmatch(size)
    if match is None:
        raise ValueError(f"Not a memory size: {size!r}")
    number, unit, binary = match.groups()
    unit = unit.upper()
    return f"{float(number) * factor:g}{unit}{'i' if binary else ''}B"


def _setting(name: str) -> str | None:
    value = os.getenv(name)
    return value.strip() if value and value.strip() else None


@dataclass(frozen=True)
class QueryBudget:
    slots: int  # queries at once
    threads: int  # DuckDB threads of the shared database (all slots)
    memory: str  # DuckDB memory limit of the shared database (all slots)


def query_budget(
    *, threads: int | None = None, memory_limit: str | None = None, cgroup: Path | None = None
) -> QueryBudget:
    """The budget from the settings, else from the container. ``threads`` and
    ``memory_limit`` (per slot) stand in for unset settings."""
    slots = int(_setting("OMNIPATH_QUERY_CONCURRENCY") or cpu_count(cgroup))
    if slots < 1:
        raise ValueError("OMNIPATH_QUERY_CONCURRENCY must be positive (query concurrency)")
    per_slot_threads = int(_setting("OMNIPATH_DUCKDB_THREADS") or threads or 1)
    if per_slot_threads < 1:
        raise ValueError("OMNIPATH_DUCKDB_THREADS must be positive")
    per_slot_memory = _setting("OMNIPATH_DUCKDB_MEMORY_LIMIT") or memory_limit
    if per_slot_memory:
        memory = scaled(per_slot_memory, slots)
    elif limit := cgroup_memory_limit(cgroup):
        memory = f"{int(limit * MEMORY_SHARE / 2**20)}MiB"
    else:
        memory = scaled(FALLBACK_MEMORY, slots)
    return QueryBudget(slots, per_slot_threads * slots, memory)
