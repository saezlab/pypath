"""COPY timing/result contract shared with explicit historical compatibility."""

from dataclasses import dataclass


@dataclass(frozen=True)
class BulkCopyResult:
    rows: int
    stage_seconds: float
    copy_seconds: float
    chunks: int
    bytes: int
