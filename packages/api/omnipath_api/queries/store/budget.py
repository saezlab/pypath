"""The query budget of an API process: how many queries run at once, and the DuckDB
threads and memory of the database they share, from the CPUs and memory its container
gets.

One query slot per CPU, and one DuckDB thread for the whole database: DuckDB runs a
query on the thread that calls it (plus the database's threads beyond the first, which
all queries share), so seven slots still run seven queries on seven cores, each without
spreading its small steps over threads that other queries need. In load tests on 7 CPUs
and 13 GB (the explorer's extended journey, 64 people): seven slots on a database of
seven threads answered the slowest 5% of steps in 0.38 s at 145 ms CPU per request, on
a database of one thread in 0.32 s at 84 ms; three slots on six threads took 1.07 s.
13 GB and 32 GB answered alike. DuckDB gets half of the container's memory, the rest is
for Python, Arrow, response buffers and DuckDB's own memory outside its limit.

Each part can be set: OMNIPATH_QUERY_CONCURRENCY (slots), OMNIPATH_DUCKDB_THREADS (the
database's threads) and OMNIPATH_DUCKDB_MEMORY_LIMIT (memory per slot). An empty setting
counts as unset."""

from __future__ import annotations

from dataclasses import dataclass
import os
from pathlib import Path
import re

from omnipath_client.budget import cgroup_memory_limit, cpu_count

# The share of the container's memory DuckDB may use. DuckDB keeps the database blocks it
# reads until it reaches its limit, so under sustained load it holds all of it, and the
# rest of the process needs room beside that: on 7 CPUs, with DuckDB at its limit (9.75 GB
# of 13 GB), the other memory (Python, Arrow, response buffers, DuckDB's allocations
# outside its limit) peaked at 2.8 GB with 64 clients on one kind of request (resolving
# 5-20 identifiers) and at 4.5 GB in the programmatic mix at 64 clients after the
# explorer journey, so 75% ran out of memory. At half (6.5 GB) the same run peaked at
# 11 GB of 13 GB and answered as fast: the cache beyond that gained nothing.
MEMORY_SHARE = 0.5
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
    threads: int  # DuckDB threads of the shared database (each query also runs on its caller's)
    memory: str  # DuckDB memory limit of the shared database (all slots)


def query_budget(
    *, threads: int | None = None, memory_limit: str | None = None, cgroup: Path | None = None
) -> QueryBudget:
    """The budget from the settings, else from the container. ``threads`` (the
    database's) and ``memory_limit`` (per slot) stand in for unset settings."""
    slots = int(_setting("OMNIPATH_QUERY_CONCURRENCY") or cpu_count(cgroup))
    if slots < 1:
        raise ValueError("OMNIPATH_QUERY_CONCURRENCY must be positive (query concurrency)")
    database_threads = int(_setting("OMNIPATH_DUCKDB_THREADS") or threads or 1)
    if database_threads < 1:
        raise ValueError("OMNIPATH_DUCKDB_THREADS must be positive")
    per_slot_memory = _setting("OMNIPATH_DUCKDB_MEMORY_LIMIT") or memory_limit
    if per_slot_memory:
        memory = scaled(per_slot_memory, slots)
    elif limit := cgroup_memory_limit(cgroup):
        memory = f"{int(limit * MEMORY_SHARE / 2**20)}MiB"
    else:
        memory = scaled(FALLBACK_MEMORY, slots)
    return QueryBudget(slots, database_threads, memory)
