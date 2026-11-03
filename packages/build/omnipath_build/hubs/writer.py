"""Streaming Parquet writer for hub mapping files."""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
import tempfile
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq

from .schema import HUB_FIELDS, HUB_SCHEMA, hub_row


class HubParquetWriter:
    """Append hub rows and stop once ``max_records`` output rows are written."""

    def __init__(
        self,
        path: str | Path,
        *,
        max_records: int | None = None,
        batch_size: int = 5_000,
        on_progress: Callable[[dict[str, Any]], None] | None = None,
        progress_stage: str | None = None,
    ) -> None:
        self.path = Path(path)
        self.max_records = max_records
        self.batch_size = batch_size
        self._on_progress = on_progress
        self._progress_stage = progress_stage
        self._buffer: list[dict[str, str]] = []
        self._count = 0
        self._writer: pq.ParquetWriter | None = None
        self.path.parent.mkdir(parents=True, exist_ok=True)
        handle = tempfile.NamedTemporaryFile(
            prefix=f".{self.path.name}.", suffix=".tmp", dir=self.path.parent, delete=False
        )
        self._staged = Path(handle.name)
        handle.close()
        self._closed = False

    @property
    def full(self) -> bool:
        return self.max_records is not None and self._count >= self.max_records

    @property
    def count(self) -> int:
        return self._count

    def add(
        self,
        source_type: str,
        source_id: str,
        hub_id: str,
        taxonomy_id: str,
        backend: str,
    ) -> bool:
        """Add one row. Returns False when the file has reached ``max_records``."""
        if self._closed:
            raise RuntimeError("Hub writer is closed")
        if self.full or not source_type or not source_id or not hub_id:
            return not self.full
        self._buffer.append(hub_row(source_type, source_id, hub_id, taxonomy_id, backend))
        self._count += 1
        if len(self._buffer) >= self.batch_size:
            self._flush()
        return not self.full

    def add_identity(
        self,
        hub_type: str,
        hub_id: str,
        taxonomy_id: str,
        backend: str,
    ) -> bool:
        return self.add(hub_type, hub_id, hub_id, taxonomy_id, backend)

    def close(self) -> int:
        """Validate and atomically publish only a successfully completed export."""
        if self._closed:
            return self._count
        try:
            self._flush()
            if self._writer is not None:
                self._writer.close()
                self._writer = None
            elif self._count == 0:
                pq.write_table(pa.Table.from_pylist([], schema=HUB_SCHEMA), self._staged)
            metadata = pq.read_metadata(self._staged)
            if metadata.num_rows != self._count or not metadata.schema.to_arrow_schema().equals(
                HUB_SCHEMA
            ):
                raise ValueError(f"Incomplete hub export: {self.path}")
            self._staged.replace(self.path)
            self._closed = True
            return self._count
        except BaseException:
            self.abort()
            raise

    def abort(self) -> None:
        """Discard partial rows, preserving the previous published hub."""
        try:
            if self._writer is not None:
                self._writer.close()
        finally:
            self._writer = None
            self._buffer.clear()
            self._staged.unlink(missing_ok=True)
            self._closed = True

    def _flush(self) -> None:
        if not self._buffer:
            return
        table = pa.Table.from_pylist(self._buffer, schema=HUB_SCHEMA).select(list(HUB_FIELDS))
        if self._writer is None:
            self._writer = pq.ParquetWriter(self._staged, HUB_SCHEMA, compression="zstd")
        self._writer.write_table(table)
        self._buffer.clear()
        if self._on_progress is not None:
            self._on_progress(
                {
                    "pipeline": "resolver",
                    "stage": self._progress_stage or "emit",
                    "status": "running",
                    "current": self._count,
                    "message": f"{self._count:,} rows",
                }
            )
