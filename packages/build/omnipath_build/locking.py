"""Process-owned build locks: stale files never prevent recovery after a crash."""

from __future__ import annotations

import fcntl
from pathlib import Path
from typing import IO


class BuildLock:
    """Nonblocking advisory lock held until close or process exit.

    The inode is deliberately retained: unlinking an advisory lock can allow two
    processes to lock different inodes for the same published resource.
    """

    def __init__(self, path: str | Path):
        self.path = Path(path)
        self.handle: IO[str] | None = None

    def __enter__(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        handle = self.path.open("a+")
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            handle.close()
            raise FileExistsError(f"Build already running: {self.path}") from exc
        self.handle = handle
        return self

    def close(self):
        if self.handle is not None:
            self.handle.close()
            self.handle = None

    def __exit__(self, *_):
        self.close()
