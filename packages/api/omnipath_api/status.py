"""Coarse serving load for the explorer's status indicator: query slots, recent response
times and the container's CPU load and memory. No host details are exposed."""

from __future__ import annotations

import os
import statistics
import time
from collections import deque
from pathlib import Path
from threading import Lock

WINDOW_SECONDS = 60
# Requests that measure the service rather than use it.
_UNTIMED = ("/health", "/status", "/api/health", "/api/status")


class RequestTimes:
    """Durations of the requests of the last minute."""

    def __init__(self, window: float = WINDOW_SECONDS):
        self.window = window
        self._samples: deque[tuple[float, float]] = deque()
        self._lock = Lock()

    def record(self, path: str, seconds: float) -> None:
        if path.endswith(_UNTIMED):
            return
        now = time.monotonic()
        with self._lock:
            self._samples.append((now, seconds))
            self._trim(now)

    def summary(self) -> tuple[int, float | None]:
        """Requests in the window and their median duration in milliseconds."""
        with self._lock:
            self._trim(time.monotonic())
            durations = [seconds for _, seconds in self._samples]
        if not durations:
            return 0, None
        return len(durations), round(statistics.median(durations) * 1000, 1)

    def _trim(self, now: float) -> None:
        while self._samples and now - self._samples[0][0] > self.window:
            self._samples.popleft()


def _read(path: str) -> str | None:
    try:
        return Path(path).read_text().strip()
    except OSError:
        return None


def memory_fraction() -> float | None:
    """The container's memory use against its limit (cgroup v2), if it has one."""
    used, limit = _read("/sys/fs/cgroup/memory.current"), _read("/sys/fs/cgroup/memory.max")
    if not used or not limit or limit == "max":
        return None
    return round(int(used) / int(limit), 3)


def cpu_load() -> float | None:
    """The one-minute load average per CPU available to the process."""
    try:
        load = os.getloadavg()[0]
    except OSError:
        return None
    cpus = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else os.cpu_count()
    return round(load / max(cpus or 1, 1), 2)


def status(pool, times: RequestTimes) -> dict:
    """``pool`` is the query pool of the API's queries (None without one)."""
    requests, median_ms = times.summary()
    running, waiting, slots = (pool.running, pool.waiting, pool.size) if pool else (0, 0, 0)
    load, memory = cpu_load(), memory_fraction()
    if waiting:
        state = "queued"
    elif (slots and running >= slots) or (load or 0) > 0.9 or (median_ms or 0) > 3000:
        state = "busy"
    else:
        state = "ok"
    return {
        "state": state,
        "queries": {"running": running, "slots": slots, "waiting": waiting},
        "requestsLastMinute": requests,
        "medianResponseMs": median_ms,
        "cpuLoad": load,
        "memory": memory,
    }
