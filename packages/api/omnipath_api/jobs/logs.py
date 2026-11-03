"""Bounded job-log reads shared by admin inspection and downloads."""

from pathlib import Path


def read_log_tail(path: Path, max_bytes: int = 2_000_000) -> str:
    with path.open("rb") as handle:
        size = handle.seek(0, 2)
        truncated = size > max_bytes
        handle.seek(max(0, size - max_bytes))
        text = handle.read(max_bytes).decode("utf-8", errors="replace")
    return ("... [truncated]\n" if truncated else "") + text
