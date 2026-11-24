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
import math
import os
from pathlib import Path
import re

CGROUP = Path("/sys/fs/cgroup")
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


def _read(path: Path) -> str | None:
    try:
        return path.read_text().strip()
    except OSError:
        return None


def _cgroups(cgroup: Path | None) -> list[Path]:
    """The directories whose limits apply: ``cgroup`` alone when given; else the process's
    own cgroup (v2, from /proc/self/cgroup: the container's root in a container, a
    systemd scope's directory outside one) and its parents, then the v1 mount root."""
    if cgroup is not None:
        return [cgroup]
    found = []
    for line in (_read(Path("/proc/self/cgroup")) or "").splitlines():
        if line.startswith("0::"):
            own = CGROUP / line[3:].strip().lstrip("/")
            found = [own, *(p for p in own.parents if p.is_relative_to(CGROUP))]
    return found or [CGROUP]


def _cpu_quota(directory: Path) -> float | None:
    if line := _read(directory / "cpu.max"):
        parts = line.split()
        if len(parts) == 2 and parts[0] != "max":
            return int(parts[0]) / int(parts[1])
        return None
    quota, period = _read(directory / "cpu" / "cpu.cfs_quota_us"), _read(directory / "cpu" / "cpu.cfs_period_us")
    if quota and period and int(quota) > 0:
        return int(quota) / int(period)
    return None


def _memory_limit(directory: Path) -> int | None:
    value = _read(directory / "memory.max") or _read(directory / "memory" / "memory.limit_in_bytes")
    if not value or value == "max":
        return None
    limit = int(value)
    return None if limit >= 2**60 else limit  # v1 writes a huge number for "no limit"


def available_cpus(cgroup: Path | None = None) -> int:
    """The CPUs the process may use: its CPU affinity, capped by the tightest CPU quota
    of its cgroups (v2 ``cpu.max``, or v1 ``cpu.cfs_quota_us``), rounded up."""
    cpus = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else os.cpu_count()
    quotas = [q for d in _cgroups(cgroup) if (q := _cpu_quota(d))]
    if quotas:
        return min(cpus or 1, max(1, math.ceil(min(quotas))))
    return cpus or 1


def memory_limit_bytes(cgroup: Path | None = None) -> int | None:
    """The tightest memory limit of the process's cgroups (v2 ``memory.max``, or v1
    ``memory.limit_in_bytes``), or None without one."""
    limits = [m for d in _cgroups(cgroup) if (m := _memory_limit(d))]
    return min(limits) if limits else None


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
    slots = int(_setting("OMNIPATH_QUERY_CONCURRENCY") or available_cpus(cgroup))
    if slots < 1:
        raise ValueError("OMNIPATH_QUERY_CONCURRENCY must be positive (query concurrency)")
    per_slot_threads = int(_setting("OMNIPATH_DUCKDB_THREADS") or threads or 1)
    if per_slot_threads < 1:
        raise ValueError("OMNIPATH_DUCKDB_THREADS must be positive")
    per_slot_memory = _setting("OMNIPATH_DUCKDB_MEMORY_LIMIT") or memory_limit
    if per_slot_memory:
        memory = scaled(per_slot_memory, slots)
    elif limit := memory_limit_bytes(cgroup):
        memory = f"{int(limit * MEMORY_SHARE / 2**20)}MiB"
    else:
        memory = scaled(FALLBACK_MEMORY, slots)
    return QueryBudget(slots, per_slot_threads * slots, memory)
