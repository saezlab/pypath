"""Bound active and idle DuckDB databases across HTTP request threads."""

from contextlib import contextmanager
from math import isfinite
from queue import Empty, LifoQueue
from threading import BoundedSemaphore, Lock


class QueryCapacityError(RuntimeError):
    """No query slot became available within the configured wait."""


class QueryPool:
    def __init__(self, factory, generation, *, size: int, timeout: float):
        if size < 1 or not isfinite(timeout) or timeout < 0:
            raise ValueError("Query concurrency must be positive and wait timeout nonnegative")
        self.factory, self.generation = factory, generation
        self.timeout = timeout
        self.size = size
        self._slots = BoundedSemaphore(size)
        self._idle = LifoQueue(size)
        # Queries holding a slot and requests waiting for one, for /status.
        self.running = 0
        self.waiting = 0
        self._lock = Lock()
        self._connections = set()
        self._closed = False

    def _discard(self, connection):
        try:
            connection.close()
        finally:
            with self._lock:
                self._connections.discard(connection)

    @contextmanager
    def acquire(self):
        with self._lock:
            self.waiting += 1
        try:
            acquired = self._slots.acquire(timeout=self.timeout)
        finally:
            with self._lock:
                self.waiting -= 1
        if not acquired:
            raise QueryCapacityError("Query capacity is busy; retry shortly")
        with self._lock:
            self.running += 1
        connection = None
        try:
            if self._closed:
                raise RuntimeError("Query pool is closed")
            generation = self.generation()
            try:
                previous_generation, connection = self._idle.get_nowait()
            except Empty:
                pass
            else:
                if previous_generation != generation:
                    self._discard(connection)
                    connection = None
            if connection is None:
                connection = self.factory()
                with self._lock:
                    self._connections.add(connection)
            yield connection
        finally:
            try:
                if connection is not None:
                    if self._closed or generation != self.generation():
                        self._discard(connection)
                    else:
                        self._idle.put_nowait((generation, connection))
            finally:
                with self._lock:
                    self.running -= 1
                self._slots.release()

    def close(self):
        """Call after draining request workers, as in the ASGI lifespan hook."""
        self._closed = True
        with self._lock:
            connections = tuple(self._connections)
        for connection in connections:
            self._discard(connection)
        while not self._idle.empty():
            self._idle.get_nowait()
