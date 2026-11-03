"""Parquet writing helpers for resolver lookup materialization."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import polars as pl
import pyarrow as pa
import pyarrow.parquet as pq


def _arrow_type(dtype: pl.DataType) -> pa.DataType:
    if dtype == pl.Utf8:
        return pa.string()
    if dtype == pl.UInt32:
        return pa.uint32()
    if dtype == pl.UInt64:
        return pa.uint64()
    raise TypeError(f'Unsupported parquet dtype conversion: {dtype!r}')


def arrow_schema(schema: dict[str, pl.DataType]) -> pa.Schema:
    """Convert a small Polars schema mapping to a PyArrow schema."""

    return pa.schema([(name, _arrow_type(dtype)) for name, dtype in schema.items()])


def write_parquet_from_dict_rows(
    rows: Iterable[dict],
    schema: dict[str, pl.DataType],
    path: str | Path,
    chunk_size: int = 100_000,
) -> int:
    """Write dictionary rows to parquet in bounded batches."""

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    pa_schema = arrow_schema(schema)
    field_names = list(schema)
    writer: pq.ParquetWriter | None = None
    count = 0
    chunk: list[dict] = []

    try:
        for row in rows:
            chunk.append({name: row.get(name) for name in field_names})
            if len(chunk) >= chunk_size:
                if writer is None:
                    writer = pq.ParquetWriter(path, pa_schema)
                writer.write_table(pa.Table.from_pylist(chunk, schema=pa_schema))
                count += len(chunk)
                chunk = []

        if chunk:
            if writer is None:
                writer = pq.ParquetWriter(path, pa_schema)
            writer.write_table(pa.Table.from_pylist(chunk, schema=pa_schema))
            count += len(chunk)
        elif writer is None:
            empty = {
                name: pa.array([], type=pa_schema.field(name).type)
                for name in field_names
            }
            pq.write_table(pa.Table.from_pydict(empty, schema=pa_schema), path)
    finally:
        if writer is not None:
            writer.close()

    return count
