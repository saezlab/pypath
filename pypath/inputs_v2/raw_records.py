"""Materialize parser-emitted raw records as content-addressed parquet."""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from datetime import UTC, datetime
import hashlib
import json
import os
from pathlib import Path
import shutil
import tempfile
import time
from typing import Any

import duckdb
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.dataset as ds
import pyarrow.parquet as pq

try:
    from omnipath_build.pipeline.progress import update_phase
except ImportError:  # pragma: no cover - pypath can run without omnipath_build
    def update_phase(detail: str = '') -> None:
        return None

PREPARSE_VERSION = 'parser_dict_columns_v1'
RAW_RECORD_BUCKET_COUNT = 4096
RAW_RECORD_PART_COUNT = 128
RAW_RECORD_MIN_PART_SIZE_BYTES = 200 * 1024 * 1024
RAW_RECORD_ID_BUCKET_STRIDE = 1_000_000_000_000
METADATA_COLUMNS = {
    '_source',
    '_dataset',
    '_raw_record_key',
    '_raw_record_id',
    'raw_record_bucket',
    'raw_record_part',
}


@dataclass(frozen=True)
class RawSnapshot:
    source: str
    dataset: str
    snapshot_id: str
    records_path: Path
    delta_path: Path
    manifest_path: Path


@dataclass(frozen=True)
class RawRecordProvenance:
    source: str
    dataset: str
    snapshot_id: str
    raw_record_key: str
    raw_record_id: int
    raw_record_bucket: int
    raw_record_part: int


@dataclass(frozen=True)
class ProvenancedRecord:
    record: Any
    provenance: RawRecordProvenance


def default_raw_records_root() -> Path:
    return Path(
        os.environ.get(
            'OMNIPATH_BRONZE_ROOT',
            os.environ.get('OMNIPATH_RAW_RECORDS_ROOT', 'data/bronze'),
        )
    )


def _log_bronze(source: str, dataset: str, message: str) -> None:
    print(f'[bronze:{source}.{dataset}] {message}', flush=True)


def materialize_raw_records(
    *,
    records: Iterable[dict[str, Any]],
    source: str,
    dataset: str,
    output_root: Path | None = None,
    batch_size: int = 50_000,
    accept: bool = False,
    download_fingerprint: dict[str, Any] | None = None,
    parser_contract: dict[str, Any] | None = None,
) -> RawSnapshot:
    """Write a complete raw-record snapshot and delta against latest."""
    output_root = output_root or default_raw_records_root()
    dataset_dir = output_root / source / dataset
    snapshot_id = datetime.now(UTC).strftime('%Y%m%dT%H%M%S%fZ')
    snapshot_dir = dataset_dir / snapshot_id
    records_path = snapshot_dir / 'records'
    delta_path = snapshot_dir / 'delta'
    manifest_path = snapshot_dir / 'manifest.json'

    latest = _read_latest(dataset_dir)
    old_snapshot_id = latest.get('snapshot_id') if latest else None
    old_records = Path(latest['records_path']) if latest else None
    if old_records is not None and not old_records.exists():
        old_records = None
        old_snapshot_id = None

    if snapshot_dir.exists():
        shutil.rmtree(snapshot_dir)
    snapshot_dir.mkdir(parents=True, exist_ok=True)

    started = datetime.now(UTC)
    phase_started = time.perf_counter()
    previous = str(old_snapshot_id) if old_snapshot_id else '-'
    _log_bronze(
        source,
        dataset,
        f'preparse start -> snapshot={snapshot_id} previous={previous} '
        f'output={snapshot_dir}',
    )
    update_phase(f'preparse snapshot={snapshot_id}')
    with tempfile.TemporaryDirectory(dir=snapshot_dir.parent) as tmpdir:
        tmp_unassigned_records = Path(tmpdir) / 'records.unassigned'
        tmp_records = Path(tmpdir) / 'records'
        tmp_delta = Path(tmpdir) / 'delta'
        _log_bronze(source, dataset, 'raw parser materialization started')
        update_phase('raw parser materialization')
        stats = _write_records(
            records,
            output_path=tmp_unassigned_records,
            source=source,
            dataset=dataset,
            snapshot_id=snapshot_id,
            batch_size=batch_size,
        )
        _log_bronze(
            source,
            dataset,
            f'raw parser materialization done -> rows={stats.get("rows", 0):,}',
        )
        update_phase(f'raw parser materialized rows={stats.get("rows", 0):,}')
        input_bytes = _parquet_size_bytes(tmp_unassigned_records)
        raw_record_part_count = _effective_part_count(
            input_bytes,
            max_part_count=RAW_RECORD_PART_COUNT,
            min_part_size_bytes=RAW_RECORD_MIN_PART_SIZE_BYTES,
        )
        _log_bronze(
            source,
            dataset,
            'raw record id assignment started -> '
            f'parts={raw_record_part_count}/{RAW_RECORD_PART_COUNT} '
            f'min_part_size={RAW_RECORD_MIN_PART_SIZE_BYTES} input_bytes={input_bytes}',
        )
        update_phase('raw record id assignment')
        id_stats = _write_records_with_ids(
            new_records=tmp_unassigned_records,
            old_records=old_records,
            output_path=tmp_records,
            part_count=raw_record_part_count,
        )
        _log_bronze(
            source,
            dataset,
            'raw record id assignment done -> '
            f'distinct_ids={id_stats.get("distinct_raw_record_ids", 0):,}',
        )
        _log_bronze(source, dataset, 'delta computation started')
        update_phase('delta computation')
        delta_stats = _write_delta(
            new_records=tmp_records,
            old_records=old_records,
            output_path=tmp_delta,
            old_snapshot_id=str(old_snapshot_id) if old_snapshot_id else None,
            new_snapshot_id=snapshot_id,
            part_count=raw_record_part_count,
        )
        _log_bronze(
            source,
            dataset,
            f'delta computation done -> {delta_stats.get("delta_keys_by_type", {})}',
        )
        shutil.move(str(tmp_records), records_path)
        shutil.move(str(tmp_delta), delta_path)

    manifest = {
        'source': source,
        'dataset': dataset,
        'snapshot_id': snapshot_id,
        'previous_snapshot_id': old_snapshot_id,
        'created_at': started.isoformat(),
        'completed_at': datetime.now(UTC).isoformat(),
        'preparse_version': PREPARSE_VERSION,
        'records_path': str(records_path),
        'delta_path': str(delta_path),
        'bucket_algorithm': 'stable_u64_sha256_mod_v1',
        'raw_record_bucket_count': RAW_RECORD_BUCKET_COUNT,
        'raw_record_part_count': raw_record_part_count,
        'requested_raw_record_part_count': RAW_RECORD_PART_COUNT,
        'raw_record_min_part_size_bytes': RAW_RECORD_MIN_PART_SIZE_BYTES,
        'raw_record_partition_input_bytes': input_bytes,
        'download_fingerprint': download_fingerprint,
        'parser_contract': parser_contract,
        **stats,
        **id_stats,
        **delta_stats,
    }
    _write_json(manifest_path, manifest)

    snapshot = RawSnapshot(
        source=source,
        dataset=dataset,
        snapshot_id=snapshot_id,
        records_path=records_path,
        delta_path=delta_path,
        manifest_path=manifest_path,
    )

    if accept:
        accept_raw_snapshot(snapshot)

    elapsed = time.perf_counter() - phase_started
    _log_bronze(
        source,
        dataset,
        f'preparse done -> rows={manifest.get("rows", 0):,} '
        f'delta={manifest.get("delta_keys_by_type", {})} '
        f'elapsed={elapsed:.1f}s',
    )
    update_phase(f'preparse done rows={manifest.get("rows", 0):,}')
    return snapshot


def reuse_raw_snapshot_if_unchanged(
    *,
    source: str,
    dataset: str,
    output_root: Path | None = None,
    download_fingerprint: dict[str, Any] | None = None,
    parser_contract: dict[str, Any] | None = None,
) -> RawSnapshot | None:
    """Create an empty-delta raw snapshot when cached download inputs are unchanged."""
    if download_fingerprint is None or parser_contract is None:
        return None

    output_root = output_root or default_raw_records_root()
    dataset_dir = output_root / source / dataset
    latest = _read_latest(dataset_dir)
    if latest is None:
        return None

    previous_manifest_path = Path(latest.get('manifest_path', ''))
    records_path = Path(latest.get('records_path', ''))
    if not previous_manifest_path.exists() or not records_path.exists():
        return None

    previous_manifest = json.loads(previous_manifest_path.read_text())
    if previous_manifest.get('download_fingerprint') != download_fingerprint:
        return None
    if previous_manifest.get('parser_contract') != parser_contract:
        return None

    snapshot_id = datetime.now(UTC).strftime('%Y%m%dT%H%M%S%fZ')
    snapshot_dir = dataset_dir / snapshot_id
    delta_path = snapshot_dir / 'delta'
    manifest_path = snapshot_dir / 'manifest.json'
    snapshot_dir.mkdir(parents=True, exist_ok=True)

    raw_record_part_count = int(previous_manifest.get('raw_record_part_count') or 1)
    _write_empty_delta(delta_path, part_count=raw_record_part_count)
    manifest = {
        'source': source,
        'dataset': dataset,
        'snapshot_id': snapshot_id,
        'previous_snapshot_id': latest.get('snapshot_id'),
        'created_at': datetime.now(UTC).isoformat(),
        'completed_at': datetime.now(UTC).isoformat(),
        'preparse_version': PREPARSE_VERSION,
        'records_path': str(records_path),
        'delta_path': str(delta_path),
        'bucket_algorithm': 'stable_u64_sha256_mod_v1',
        'raw_record_bucket_count': RAW_RECORD_BUCKET_COUNT,
        'raw_record_part_count': raw_record_part_count,
        'requested_raw_record_part_count': RAW_RECORD_PART_COUNT,
        'raw_record_min_part_size_bytes': RAW_RECORD_MIN_PART_SIZE_BYTES,
        'download_fingerprint': download_fingerprint,
        'parser_contract': parser_contract,
        'reused_from_snapshot_id': latest.get('snapshot_id'),
        'reused_without_preparse': True,
        'delta_keys_by_type': {},
        'rows': int(previous_manifest.get('rows', 0) or 0),
        'min_raw_record_id': int(previous_manifest.get('min_raw_record_id', 0) or 0),
        'max_raw_record_id': int(previous_manifest.get('max_raw_record_id', 0) or 0),
        'distinct_raw_record_ids': int(previous_manifest.get('distinct_raw_record_ids', 0) or 0),
    }
    _write_json(manifest_path, manifest)
    _log_bronze(
        source,
        dataset,
        f'download unchanged; reused raw snapshot {latest.get("snapshot_id")} '
        f'-> snapshot={snapshot_id}',
    )
    update_phase(f'preparse skipped snapshot={snapshot_id}')
    return RawSnapshot(
        source=source,
        dataset=dataset,
        snapshot_id=snapshot_id,
        records_path=records_path,
        delta_path=delta_path,
        manifest_path=manifest_path,
    )


def accept_raw_snapshot(snapshot: RawSnapshot) -> None:
    """Advance latest.json to a successfully processed raw snapshot.

    Full records are materialized under the snapshot directory while downstream
    silver processing runs. On accept, move that complete records file to a
    single mutable state location so unchanged raw payload rows are not kept in
    every snapshot directory. Immutable snapshot directories retain the compact
    delta and manifest.
    """
    dataset_dir = snapshot.records_path.parent.parent
    state_dir = dataset_dir / 'state'
    state_records_path = state_dir / 'records'
    state_dir.mkdir(parents=True, exist_ok=True)

    accepted_records_path = snapshot.records_path
    if snapshot.records_path != state_records_path:
        tmp_state_records_path = state_dir / 'records.tmp'
        if tmp_state_records_path.exists():
            shutil.rmtree(tmp_state_records_path)
        shutil.move(str(snapshot.records_path), tmp_state_records_path)
        if state_records_path.exists():
            shutil.rmtree(state_records_path)
        tmp_state_records_path.replace(state_records_path)
        accepted_records_path = state_records_path

    _rewrite_manifest_records_path(snapshot.manifest_path, accepted_records_path)
    _write_json(
        dataset_dir / 'latest.json',
        {
            'source': snapshot.source,
            'dataset': snapshot.dataset,
            'snapshot_id': snapshot.snapshot_id,
            'records_path': str(accepted_records_path),
            'delta_path': str(snapshot.delta_path),
            'manifest_path': str(snapshot.manifest_path),
            'updated_at': datetime.now(UTC).isoformat(),
        },
    )


def iter_raw_record_dicts(
    records_path: Path,
    *,
    keys: set[str] | None = None,
    include_metadata: bool = False,
) -> Iterator[dict[str, Any]]:
    """Yield mapper-input dictionaries from a partitioned raw records dataset."""
    key_filter = keys
    dataset = ds.dataset(records_path, format='parquet')
    for batch in dataset.to_batches(batch_size=10_000):
        table = pa.Table.from_batches([batch], schema=batch.schema)
        for row in table.to_pylist():
            if key_filter is not None and row.get('_raw_record_key') not in key_filter:
                continue
            if include_metadata:
                yield row
            else:
                yield {k: v for k, v in row.items() if k not in METADATA_COLUMNS}


def iter_changed_raw_record_dicts(
    records_path: Path,
    delta_path: Path,
    *,
    include_metadata: bool = False,
) -> Iterator[dict[str, Any]]:
    """Yield changed current records by joining partitioned records to delta."""
    con = duckdb.connect()
    try:
        reader = con.execute(
            f"""
            SELECT r.*
            FROM {_read_records_sql(records_path)} r
            JOIN {_read_records_sql(delta_path)} d USING (_raw_record_key)
            WHERE d._change_type = 'added'
            ORDER BY r.raw_record_bucket, r._raw_record_key
            """
        ).fetch_record_batch(rows_per_batch=10_000)
        while True:
            try:
                batch = reader.read_next_batch()
            except StopIteration:
                break
            table = pa.Table.from_batches([batch])
            for row in table.to_pylist():
                if include_metadata:
                    yield row
                else:
                    yield {k: v for k, v in row.items() if k not in METADATA_COLUMNS}
    finally:
        con.close()


def changed_keys(delta_path: Path) -> Iterator[str]:
    if not delta_path.exists():
        return
    dataset = ds.dataset(delta_path, format='parquet')
    for batch in dataset.to_batches(columns=['_raw_record_key'], batch_size=10_000):
        for value in batch.column('_raw_record_key').to_pylist():
            if value is None:
                continue
            yield str(value)


def changed_key_count(delta_path: Path) -> int:
    if not delta_path.exists():
        return 0
    con = duckdb.connect()
    try:
        return int(con.execute(f"""
            SELECT count(DISTINCT _raw_record_key)
            FROM {_read_records_sql(delta_path)}
            WHERE _raw_record_key IS NOT NULL
        """).fetchone()[0] or 0)
    finally:
        con.close()


def changed_key_counts_by_type(delta_path: Path) -> dict[str, int]:
    if not delta_path.exists():
        return {'added': 0, 'removed': 0}
    con = duckdb.connect()
    try:
        rows = con.execute(f"""
            SELECT _change_type, count(DISTINCT _raw_record_key) AS count
            FROM {_read_records_sql(delta_path)}
            WHERE _raw_record_key IS NOT NULL
            GROUP BY _change_type
        """).fetchall()
    finally:
        con.close()
    counts = {str(change_type): int(count) for change_type, count in rows}
    return {
        'added': counts.get('added', 0),
        'removed': counts.get('removed', 0),
    }


def _write_empty_delta(output_path: Path, *, part_count: int = 1) -> None:
    if output_path.exists():
        shutil.rmtree(output_path)
    output_path.mkdir(parents=True, exist_ok=True)
    table = pa.Table.from_arrays(
        [
            pa.array([], type=pa.string()),
            pa.array([], type=pa.int64()),
            pa.array([], type=pa.int64()),
            pa.array([], type=pa.int64()),
            pa.array([], type=pa.string()),
        ],
        names=[
            '_raw_record_key',
            '_raw_record_id',
            'raw_record_bucket',
            'raw_record_part',
            '_change_type',
        ],
    )
    for part in range(max(1, part_count)):
        part_dir = output_path / f'part={part:05d}'
        part_dir.mkdir(parents=True, exist_ok=True)
        pq.write_table(table, part_dir / 'data.parquet', compression='zstd')


def _write_records(
    records: Iterable[dict[str, Any]],
    *,
    output_path: Path,
    source: str,
    dataset: str,
    snapshot_id: str,
    batch_size: int,
) -> dict[str, Any]:
    if output_path.exists():
        shutil.rmtree(output_path)
    output_path.mkdir(parents=True, exist_ok=True)
    writers: dict[int, pq.ParquetWriter] = {}
    batch: list[dict[str, Any]] = []
    schema_names: list[str] | None = None
    schema: pa.Schema | None = None
    rows = 0
    last_logged_rows = 0
    started = time.perf_counter()

    def flush() -> None:
        nonlocal batch, schema_names, schema, last_logged_rows
        if not batch:
            return
        if schema_names is None:
            seen: dict[str, None] = {}
            for row in batch:
                for name in row:
                    seen.setdefault(name, None)
            schema_names = list(seen)
        normalized = [
            {name: _stringify_if_unsupported(row.get(name)) for name in schema_names}
            for row in batch
        ]
        table = pa.Table.from_pylist(normalized)
        if schema is None:
            schema = table.schema
        else:
            table = table.cast(schema, safe=False)
        for part in sorted(set(table.column('raw_record_part').to_pylist())):
            if part is None:
                continue
            part_int = int(part)
            part_table = table.filter(
                pc.equal(table.column('raw_record_part'), pa.scalar(part_int, type=pa.int64()))
            )
            writer = writers.get(part_int)
            if writer is None:
                part_dir = output_path / f'part={part_int:05d}'
                part_dir.mkdir(parents=True, exist_ok=True)
                writer = pq.ParquetWriter(
                    part_dir / 'data.parquet',
                    schema,
                    compression='zstd',
                    use_dictionary=True,
                )
                writers[part_int] = writer
            writer.write_table(part_table)
        batch = []
        if rows - last_logged_rows >= batch_size:
            elapsed = time.perf_counter() - started
            _log_bronze(
                source,
                dataset,
                f'raw parser materialized {rows:,} rows in {elapsed:.1f}s',
            )
            update_phase(f'raw parser materialized rows={rows:,}')
            last_logged_rows = rows

    try:
        for record in records:
            clean = _clean_record(record)
            key = _record_hash(clean)
            bucket = _stable_bucket(key, RAW_RECORD_BUCKET_COUNT)
            part = _stable_part_from_bucket(bucket, RAW_RECORD_BUCKET_COUNT, RAW_RECORD_PART_COUNT)
            row = {
                '_source': source,
                '_dataset': dataset,
                '_raw_record_key': key,
                'raw_record_bucket': bucket,
                'raw_record_part': part,
                **clean,
            }
            if schema_names is not None:
                missing = set(row) - set(schema_names)
                if missing:
                    raise ValueError(
                        'Raw parser emitted new columns after first batch: '
                        f'{sorted(missing)}. Use a stable parser schema.'
                    )
            batch.append(row)
            rows += 1
            if len(batch) >= batch_size:
                flush()
        flush()
    finally:
        if schema is None:
            schema = pa.schema([
                pa.field('_source', pa.string()),
                pa.field('_dataset', pa.string()),
                pa.field('_raw_record_key', pa.string()),
                pa.field('raw_record_bucket', pa.int64()),
                pa.field('raw_record_part', pa.int64()),
            ])
        empty = pa.Table.from_pylist([], schema=schema)
        if not writers:
            part_dir = output_path / 'part=00000'
            part_dir.mkdir(parents=True, exist_ok=True)
            pq.write_table(empty, part_dir / 'data.parquet', compression='zstd')
        for writer in writers.values():
            writer.close()

    duplicate_stats = _duplicate_stats(output_path) if rows else {}
    return {'rows': rows, **duplicate_stats}


def _clean_record(record: dict[str, Any]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for key, value in record.items():
        name = str(key) if key is not None else 'column'
        if name in METADATA_COLUMNS:
            name = f'raw{name}'
        out[name] = value
    return out


def _record_hash(record: dict[str, Any]) -> str:
    digest = hashlib.blake2b(digest_size=32)
    for key in sorted(record):
        digest.update(key.encode('utf-8', errors='surrogatepass'))
        digest.update(b'\x1f')
        digest.update(_canonical_bytes(record[key]))
        digest.update(b'\x1e')
    return digest.hexdigest()


def _canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        _canonical_value(value),
        sort_keys=True,
        separators=(',', ':'),
        ensure_ascii=False,
    ).encode('utf-8', errors='surrogatepass')


def _canonical_value(value: Any) -> Any:
    if value is None:
        return None
    if isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, bytes):
        return value.decode('utf-8', errors='replace')
    if isinstance(value, (list, tuple)):
        return [_canonical_value(item) for item in value]
    if isinstance(value, dict):
        return {str(k): _canonical_value(v) for k, v in sorted(value.items())}
    return str(value)


def _stringify_if_unsupported(value: Any) -> Any:
    # Keep common parquet-native values. Fall back to deterministic JSON for
    # exotic Python objects so every raw parser can be materialized.
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, bytes):
        return value.decode('utf-8', errors='replace')
    if isinstance(value, list):
        return [_stringify_if_unsupported(item) for item in value]
    if isinstance(value, tuple):
        return [_stringify_if_unsupported(item) for item in value]
    if isinstance(value, dict):
        return {str(k): _stringify_if_unsupported(v) for k, v in value.items()}
    return str(value)


def _write_records_with_ids(
    *,
    new_records: Path,
    old_records: Path | None,
    output_path: Path,
    part_count: int,
) -> dict[str, Any]:
    """Copy new records while reusing/assigning bucket-local raw-record IDs."""
    if output_path.exists():
        shutil.rmtree(output_path)
    output_path.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect()
    try:
        part_expr = _raw_record_part_sql(part_count)
        part_expr_n = _raw_record_part_sql(part_count, table_alias='n')
        for part in range(part_count):
            part_dir = output_path / f'part={part:05d}'
            part_dir.mkdir(parents=True, exist_ok=True)
            output_sql = _sql_literal(str(part_dir / 'data.parquet'))
            old_map_sql = _old_raw_record_id_map_sql(old_records, part, part_count)
            con.execute(
                f"""
                COPY (
                    WITH
                    old_map AS ({old_map_sql}),
                    new_keys AS (
                        SELECT DISTINCT
                            _raw_record_key,
                            raw_record_bucket,
                            {part_expr} AS raw_record_part
                        FROM {_read_records_sql(new_records)}
                        WHERE {part_expr} = {part}
                    ),
                    max_old AS (
                        SELECT
                            raw_record_bucket,
                            max(_raw_record_id % {RAW_RECORD_ID_BUCKET_STRIDE}) AS max_local_id
                        FROM old_map
                        GROUP BY raw_record_bucket
                    ),
                    added_map AS (
                        SELECT
                            k._raw_record_key,
                            (
                                k.raw_record_bucket * {RAW_RECORD_ID_BUCKET_STRIDE}
                                + coalesce(m.max_local_id, 0)
                                + row_number() OVER (
                                    PARTITION BY k.raw_record_bucket
                                    ORDER BY k._raw_record_key
                                )
                            )::BIGINT AS _raw_record_id,
                            k.raw_record_bucket,
                            k.raw_record_part
                        FROM new_keys k
                        LEFT JOIN max_old m USING (raw_record_bucket)
                        WHERE k._raw_record_key NOT IN (
                            SELECT _raw_record_key FROM old_map
                        )
                    ),
                    id_map AS (
                        SELECT
                            _raw_record_key,
                            _raw_record_id,
                            raw_record_bucket,
                            raw_record_part
                        FROM old_map
                        WHERE _raw_record_key IN (
                            SELECT _raw_record_key FROM new_keys
                        )
                        UNION ALL
                        SELECT
                            _raw_record_key,
                            _raw_record_id,
                            raw_record_bucket,
                            raw_record_part
                        FROM added_map
                    )
                    SELECT
                        n._source,
                        n._dataset,
                        n._raw_record_key,
                        CAST(m._raw_record_id AS BIGINT) AS _raw_record_id,
                        n.raw_record_bucket,
                        {part_expr_n} AS raw_record_part,
                        n.* EXCLUDE (
                            _source,
                            _dataset,
                            _raw_record_key,
                            raw_record_bucket,
                            raw_record_part
                        )
                    FROM {_read_records_sql(new_records)} AS n
                    JOIN id_map AS m USING (_raw_record_key)
                    WHERE {part_expr_n} = {part}
                    ORDER BY n.raw_record_bucket, n._raw_record_key
                ) TO {output_sql} (FORMAT PARQUET, COMPRESSION ZSTD)
                """
            )
        min_id, max_id, distinct_ids = con.execute(
            f"""
            SELECT min(_raw_record_id), max(_raw_record_id), count(DISTINCT _raw_record_id)
            FROM {_read_records_sql(output_path)}
            """
        ).fetchone()
    finally:
        con.close()
    return {
        'min_raw_record_id': int(min_id or 0),
        'max_raw_record_id': int(max_id or 0),
        'distinct_raw_record_ids': int(distinct_ids or 0),
    }


def _write_delta(
    *,
    new_records: Path,
    old_records: Path | None,
    output_path: Path,
    old_snapshot_id: str | None,
    new_snapshot_id: str,
    part_count: int,
) -> dict[str, Any]:
    if output_path.exists():
        shutil.rmtree(output_path)
    output_path.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect()
    try:
        part_expr = _raw_record_part_sql(part_count)
        for part in range(part_count):
            part_dir = output_path / f'part={part:05d}'
            part_dir.mkdir(parents=True, exist_ok=True)
            output_sql = _sql_literal(str(part_dir / 'data.parquet'))
            old_map_sql = _old_raw_record_id_map_sql(old_records, part, part_count)
            con.execute(
                f"""
                COPY (
                    WITH
                    old_keys AS ({old_map_sql}),
                    new_keys AS (
                        SELECT DISTINCT
                            _raw_record_key,
                            _raw_record_id,
                            raw_record_bucket,
                            {part_expr} AS raw_record_part
                        FROM {_read_records_sql(new_records)}
                        WHERE {part_expr} = {part}
                    )
                    SELECT
                        _raw_record_key,
                        _raw_record_id,
                        raw_record_bucket,
                        raw_record_part,
                        'added' AS _change_type
                    FROM new_keys
                    WHERE _raw_record_key NOT IN (SELECT _raw_record_key FROM old_keys)
                    UNION ALL
                    SELECT
                        _raw_record_key,
                        _raw_record_id,
                        raw_record_bucket,
                        raw_record_part,
                        'removed' AS _change_type
                    FROM old_keys
                    WHERE _raw_record_key NOT IN (SELECT _raw_record_key FROM new_keys)
                    ORDER BY raw_record_bucket, _raw_record_key, _change_type
                ) TO {output_sql} (FORMAT PARQUET, COMPRESSION ZSTD)
                """
            )
        stats = con.execute(
            f"""
            SELECT _change_type, count(*) AS keys
            FROM {_read_records_sql(output_path)}
            GROUP BY _change_type
            ORDER BY _change_type
            """
        ).fetchall()
    finally:
        con.close()
    return {'delta_keys_by_type': dict(stats)}


def _old_raw_record_id_map_sql(old_records: Path | None, part: int, part_count: int) -> str:
    if old_records is None:
        return """
            SELECT
                '' AS _raw_record_key,
                CAST(NULL AS BIGINT) AS _raw_record_id,
                CAST(NULL AS BIGINT) AS raw_record_bucket,
                CAST(NULL AS BIGINT) AS raw_record_part
            WHERE false
        """
    return f"""
        SELECT
            _raw_record_key,
            min(_raw_record_id)::BIGINT AS _raw_record_id,
            min(raw_record_bucket)::BIGINT AS raw_record_bucket,
            min({_raw_record_part_sql(part_count)})::BIGINT AS raw_record_part
        FROM {_read_records_sql(old_records)}
        WHERE {_raw_record_part_sql(part_count)} = {part}
        GROUP BY _raw_record_key
    """


def _duplicate_stats(records_path: Path) -> dict[str, Any]:
    con = duckdb.connect()
    try:
        duplicate_key_count, duplicate_row_count = con.execute(
            f"""
            WITH counts AS (
                SELECT _raw_record_key, count(*) AS n
                FROM {_read_records_sql(records_path)}
                GROUP BY _raw_record_key
            )
            SELECT
                count(*) FILTER (WHERE n > 1) AS duplicate_key_count,
                coalesce(sum(n - 1) FILTER (WHERE n > 1), 0) AS duplicate_row_count
            FROM counts
            """
        ).fetchone()
    finally:
        con.close()
    return {
        'duplicate_key_count': int(duplicate_key_count or 0),
        'duplicate_row_count': int(duplicate_row_count or 0),
    }


def _rewrite_manifest_records_path(manifest_path: Path, records_path: Path) -> None:
    if not manifest_path.exists():
        return
    manifest = json.loads(manifest_path.read_text())
    manifest['records_path'] = str(records_path)
    manifest['accepted_records_path'] = str(records_path)
    manifest['accepted_at'] = datetime.now(UTC).isoformat()
    _write_json(manifest_path, manifest)


def _read_latest(dataset_dir: Path) -> dict[str, Any] | None:
    path = dataset_dir / 'latest.json'
    if not path.exists():
        return None
    return json.loads(path.read_text())


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + '\n')


def _sql_literal(value: str) -> str:
    return "'" + value.replace("'", "''") + "'"


def _read_records_sql(path: Path) -> str:
    return (
        "read_parquet("
        f"{_sql_literal(str(path / '**' / '*.parquet'))}, "
        "union_by_name=true, hive_partitioning=true)"
    )


def _raw_record_part_sql(part_count: int, *, table_alias: str | None = None) -> str:
    bucket_column = (
        f'{table_alias}.raw_record_bucket'
        if table_alias is not None
        else 'raw_record_bucket'
    )
    return (
        f"floor({bucket_column} * {part_count} / "
        f"{RAW_RECORD_BUCKET_COUNT})::BIGINT"
    )


def _parquet_size_bytes(path: Path) -> int:
    if not path.exists():
        return 0
    if path.is_file():
        return path.stat().st_size if path.suffix == '.parquet' else 0
    return sum(file.stat().st_size for file in path.rglob('*.parquet') if file.is_file())


def _effective_part_count(
    input_bytes: int,
    *,
    max_part_count: int,
    min_part_size_bytes: int,
) -> int:
    if min_part_size_bytes <= 0 or input_bytes <= 0:
        return 1
    return max(1, min(max_part_count, input_bytes // min_part_size_bytes))


def _stable_bucket(value: str, bucket_count: int) -> int:
    digest = hashlib.sha256(value.encode('utf-8')).digest()
    return int.from_bytes(digest[:8], 'big', signed=False) % bucket_count


def _stable_part_from_bucket(bucket: int, bucket_count: int, part_count: int) -> int:
    return int(bucket * part_count // bucket_count)
