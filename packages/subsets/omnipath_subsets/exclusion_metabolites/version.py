"""Version of the seed and override files, stamped on every row.

A committed file is identified by the commit that last changed the set. A file
with edits that no commit holds, or one that is not under git at all, is
identified by a hash of the content, and the report says so. The same inputs
always give the same value.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
from pathlib import Path
import subprocess


@dataclass(frozen=True)
class ListVersion:
    value: str
    uncommitted: bool


def content_hash(paths) -> str:
    sha = hashlib.sha256()
    for path in sorted(paths, key=lambda p: p.name):
        sha.update(path.name.encode() + b"\0" + path.read_bytes() + b"\0")
    return sha.hexdigest()


def _git(directory: Path, *args: str) -> str | None:
    try:
        result = subprocess.run(
            ["git", *args], cwd=directory, capture_output=True, text=True, timeout=30
        )
    except (OSError, subprocess.SubprocessError):
        return None
    return result.stdout.strip() if result.returncode == 0 else None


def list_version(paths) -> ListVersion:
    """The commit of ``paths`` if all are committed and unchanged, else a hash."""
    paths = [Path(path) for path in paths]
    directory = paths[0].parent
    names = [path.name for path in paths]
    dirty = _git(directory, "status", "--porcelain", "--", *names)
    head = _git(directory, "log", "-1", "--format=%H", "--", *names)
    if dirty == "" and head:
        return ListVersion(head, False)
    return ListVersion("uncommitted:" + content_hash(paths), True)
