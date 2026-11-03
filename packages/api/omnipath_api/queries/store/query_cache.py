"""Bounded result cache that shares simultaneous identical read queries."""

from collections import OrderedDict
from concurrent.futures import Future
from copy import deepcopy
import json
import threading


def normalized(value):
    """Facet arrays are sets; omitted and empty fields have the same meaning."""
    if isinstance(value, dict):
        return {
            k: normalized(v)
            for k, v in sorted(value.items())
            if v is not None and v != [] and v != {}
        }
    if isinstance(value, (list, tuple)):
        return sorted({json.dumps(normalized(v), sort_keys=True) for v in value})
    return value


class QueryCache:
    def __init__(self, max_entries=128, max_bytes=32 * 1024**2):
        self.max_entries, self.max_bytes = max_entries, max_bytes
        self.lock = threading.Lock()
        self.entries = OrderedDict()
        self.pending = {}
        self.bytes = 0

    def clear(self):
        with self.lock:
            self.entries.clear()
            self.bytes = 0

    def get(self, key, compute):
        with self.lock:
            if key in self.entries:
                value, _ = self.entries[key]
                self.entries.move_to_end(key)
                return deepcopy(value)
            future = self.pending.get(key)
            owner = future is None
            if owner:
                future = self.pending[key] = Future()
        if not owner:
            return deepcopy(future.result())
        try:
            value = compute()
            size = len(json.dumps(value, default=str).encode())
            with self.lock:
                if size <= self.max_bytes:
                    self.entries[key] = (deepcopy(value), size)
                    self.bytes += size
                    while len(self.entries) > self.max_entries or self.bytes > self.max_bytes:
                        _, (_, removed) = self.entries.popitem(last=False)
                        self.bytes -= removed
            future.set_result(value)
            return deepcopy(value)
        except BaseException as exc:
            future.set_exception(exc)
            raise
        finally:
            with self.lock:
                self.pending.pop(key, None)
