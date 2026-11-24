"""The defaults of the database build: DuckDB threads and memory from the machine.

The build's memory need grows with the threads, and it holds steady at about 5 GB at 3
threads for a full release (173 resources; the hosted API builds it that way in about
22 minutes), while a build with every core of a large machine and 2 GB ran out of
memory. So the defaults are: at most 3 threads (fewer CPUs, fewer threads), and 60% of
the memory the process can use, never below 2 GB. That memory is the machine's
available RAM, or the tighter limit of its container (cgroup), if any. DuckDB takes
only what it needs from the limit and writes the rest to disk.
"""

from __future__ import annotations

import math
import os
from pathlib import Path

MAX_THREADS = 3
MEMORY_SHARE = 0.6
MIN_MEMORY = 2 * 2**30

CGROUP = Path("/sys/fs/cgroup")


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


def cpu_count(cgroup: Path | None = None) -> int:
    """The CPUs the process may use: its CPU affinity, capped by the tightest CPU quota
    of its cgroups (v2 ``cpu.max``, or v1 ``cpu.cfs_quota_us``), rounded up. ``cgroup``
    names the one directory to read (for tests)."""
    cpus = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else os.cpu_count()
    quotas = [q for d in _cgroups(cgroup) if (q := _cpu_quota(d))]
    if quotas:
        return min(cpus or 1, max(1, math.ceil(min(quotas))))
    return cpus or 1


def cgroup_memory_limit(cgroup: Path | None = None) -> int | None:
    """The tightest memory limit of the process's cgroups (v2 ``memory.max``, or v1
    ``memory.limit_in_bytes``), or None without one."""
    limits = [m for d in _cgroups(cgroup) if (m := _memory_limit(d))]
    return min(limits) if limits else None


def _system_memory() -> int | None:
    """The memory available to start with: MemAvailable on Linux, else all physical RAM."""
    for line in (_read(Path("/proc/meminfo")) or "").splitlines():
        if line.startswith("MemAvailable:"):
            return int(line.split()[1]) * 1024
    try:
        return os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
    except (AttributeError, ValueError, OSError):
        return None


def available_memory() -> int | None:
    """Bytes of memory the process can use: the machine's available RAM, or the cgroup
    limit if lower; None when unknown."""
    known = [m for m in (_system_memory(), cgroup_memory_limit()) if m]
    return min(known) if known else None


def build_defaults(memory_limit: str | None = None, threads: int | None = None) -> tuple[str, int]:
    """(DuckDB memory limit, threads) for a database build; what is given stays."""
    if threads is not None and threads < 1:
        raise ValueError("threads must be positive")
    if threads is None:
        threads = min(cpu_count(), MAX_THREADS)
    if memory_limit is None:
        memory = available_memory()
        size = max(int(memory * MEMORY_SHARE), MIN_MEMORY) if memory else MIN_MEMORY * 2
        memory_limit = f"{size // 2**20}MiB"
    return memory_limit, threads


def out_of_memory_hint(memory_limit: str, threads: int) -> str:
    """What to do when a build ran out of memory."""
    return (
        f"The database build ran out of memory (limit {memory_limit}, {threads} threads). "
        "Give it more (memory_limit / --memory-limit, e.g. 6GB; a full release takes about "
        "5 GB at 3 threads), use fewer threads (threads / --threads), or install fewer "
        "resources."
    )
