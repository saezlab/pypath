"""Materialize parser-emitted raw records as content-addressed parquet."""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
import hashlib
import json
import os
import shutil
import tempfile
from typing import Any

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq

PREPARSE_VERSION = 'parser_dict_columns_v1'
METADATA_COLUMNS = {
    '_source',
    '_dataset',
    '_raw_record_key',
    '_raw_record_id',
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


def materialize_raw_records(
    *,
    records: Iterable[dict[str, Any]],
    source: str,
    dataset: str,
    output_root: Path | None = None,
    batch_size: int = 50_000,
    accept: bool = False,
) -> RawSnapshot:
    """Write a complete raw-record snapshot and delta against latest."""
    output_root = output_root or default_raw_records_root()
    dataset_dir = output_root / source / dataset
    snapshot_id = datetime.now(UTC).strftime('%Y%m%dT%H%M%S%fZ')
    snapshot_dir = dataset_dir / snapshot_id
    records_path = snapshot_dir / 'records.parquet'
    delta_path = snapshot_dir / 'delta.parquet'
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
    with tempfile.TemporaryDirectory(dir=snapshot_dir.parent) as tmpdir:
        tmp_unassigned_records = Path(tmpdir) / 'records.unassigned.parquet'
        tmp_records = Path(tmpdir) / 'records.parquet'
        tmp_delta = Path(tmpdir) / 'delta.parquet'
        stats = _write_records(
            records,
            output_path=tmp_unassigned_records,
            source=source,
            dataset=dataset,
            snapshot_id=snapshot_id,
            batch_size=batch_size,
        )
        id_stats = _write_records_with_ids(
            new_records=tmp_unassigned_records,
            old_records=old_records,
            output_path=tmp_records,
        )
        delta_stats = _write_delta(
            new_records=tmp_records,
            old_records=old_records,
            output_path=tmp_delta,
            old_snapshot_id=str(old_snapshot_id) if old_snapshot_id else None,
            new_snapshot_id=snapshot_id,
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

    return snapshot


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
    state_records_path = state_dir / 'records.parquet'
    state_dir.mkdir(parents=True, exist_ok=True)

    accepted_records_path = snapshot.records_path
    if snapshot.records_path != state_records_path:
        tmp_state_records_path = state_dir / 'records.parquet.tmp'
        if tmp_state_records_path.exists():
            tmp_state_records_path.unlink()
        shutil.move(str(snapshot.records_path), tmp_state_records_path)
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
    """Yield mapper-input dictionaries from a raw records parquet file."""
    key_filter = keys
    parquet = pq.ParquetFile(records_path)
    for batch in parquet.iter_batches(batch_size=10_000):
        table = pa.Table.from_batches([batch])
        for row in table.to_pylist():
            if key_filter is not None and row.get('_raw_record_key') not in key_filter:
                continue
            if include_metadata:
                yield row
            else:
                yield {k: v for k, v in row.items() if k not in METADATA_COLUMNS}


def changed_keys(delta_path: Path) -> set[str]:
    if not delta_path.exists():
        return set()
    table = pq.read_table(delta_path, columns=['_raw_record_key'])
    return set(table.column('_raw_record_key').to_pylist())


def _write_records(
    records: Iterable[dict[str, Any]],
    *,
    output_path: Path,
    source: str,
    dataset: str,
    snapshot_id: str,
    batch_size: int,
) -> dict[str, Any]:
    writer: pq.ParquetWriter | None = None
    batch: list[dict[str, Any]] = []
    schema_names: list[str] | None = None
    rows = 0

    def flush() -> None:
        nonlocal writer, batch, schema_names
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
        if writer is None:
            writer = pq.ParquetWriter(
                output_path,
                table.schema,
                compression='zstd',
                use_dictionary=True,
            )
        else:
            table = table.cast(writer.schema, safe=False)
        writer.write_table(table)
        batch = []

    try:
        for record in records:
            clean = _clean_record(record)
            key = _record_hash(clean)
            row = {
                '_source': source,
                '_dataset': dataset,
                '_raw_record_key': key,
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
        if writer is not None:
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
) -> dict[str, Any]:
    """Copy new records while reusing/assigning compact raw-record IDs."""
    con = duckdb.connect()
    output_sql = _sql_literal(str(output_path))
    old_map_sql = _old_raw_record_id_map_sql(old_records)
    try:
        con.execute(
            f"""
            COPY (
                WITH
                old_map AS ({old_map_sql}),
                new_keys AS (
                    SELECT DISTINCT _raw_record_key
                    FROM read_parquet(?)
                ),
                max_old AS (
                    SELECT coalesce(max(_raw_record_id), 0) AS max_id FROM old_map
                ),
                added_map AS (
                    SELECT
                        _raw_record_key,
                        (SELECT max_id FROM max_old)
                        + row_number() OVER (ORDER BY _raw_record_key) AS _raw_record_id
                    FROM new_keys
                    WHERE _raw_record_key NOT IN (SELECT _raw_record_key FROM old_map)
                ),
                id_map AS (
                    SELECT _raw_record_key, _raw_record_id FROM old_map
                    WHERE _raw_record_key IN (SELECT _raw_record_key FROM new_keys)
                    UNION ALL
                    SELECT _raw_record_key, _raw_record_id FROM added_map
                )
                SELECT
                    n._source,
                    n._dataset,
                    n._raw_record_key,
                    CAST(m._raw_record_id AS BIGINT) AS _raw_record_id,
                    n.* EXCLUDE (_source, _dataset, _raw_record_key)
                FROM read_parquet(?) AS n
                JOIN id_map AS m USING (_raw_record_key)
            ) TO {output_sql} (FORMAT PARQUET, COMPRESSION ZSTD)
            """,
            [str(new_records), str(new_records)],
        )
        min_id, max_id, distinct_ids = con.execute(
            """
            SELECT min(_raw_record_id), max(_raw_record_id), count(DISTINCT _raw_record_id)
            FROM read_parquet(?)
            """,
            [str(output_path)],
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
) -> dict[str, Any]:
    con = duckdb.connect()
    output_sql = _sql_literal(str(output_path))
    old_map_sql = _old_raw_record_id_map_sql(old_records)
    try:
        if old_records is None:
            con.execute(
                f"""
                COPY (
                    SELECT DISTINCT
                        _raw_record_key,
                        _raw_record_id,
                        'added' AS _change_type
                    FROM read_parquet(?)
                ) TO {output_sql} (FORMAT PARQUET, COMPRESSION ZSTD)
                """,
                [str(new_records)],
            )
        else:
            con.execute(
                f"""
                COPY (
                    WITH
                    old_keys AS ({old_map_sql}),
                    new_keys AS (
                        SELECT DISTINCT _raw_record_key, _raw_record_id
                        FROM read_parquet(?)
                    )
                    SELECT _raw_record_key, _raw_record_id, 'added' AS _change_type
                    FROM new_keys
                    WHERE _raw_record_key NOT IN (SELECT _raw_record_key FROM old_keys)
                    UNION ALL
                    SELECT _raw_record_key, _raw_record_id, 'removed' AS _change_type
                    FROM old_keys
                    WHERE _raw_record_key NOT IN (SELECT _raw_record_key FROM new_keys)
                ) TO {output_sql} (FORMAT PARQUET, COMPRESSION ZSTD)
                """,
                [str(new_records)],
            )
        stats = con.execute(
            """
            SELECT _change_type, count(*) AS keys
            FROM read_parquet(?)
            GROUP BY _change_type
            ORDER BY _change_type
            """,
            [str(output_path)],
        ).fetchall()
    finally:
        con.close()
    return {'delta_keys_by_type': dict(stats)}


def _old_raw_record_id_map_sql(old_records: Path | None) -> str:
    if old_records is None:
        return "SELECT '' AS _raw_record_key, CAST(NULL AS BIGINT) AS _raw_record_id WHERE false"

    old_path = _sql_literal(str(old_records))
    schema = pq.read_schema(old_records)
    if '_raw_record_id' in schema.names:
        return f"""
            SELECT _raw_record_key, min(_raw_record_id)::BIGINT AS _raw_record_id
            FROM read_parquet({old_path})
            GROUP BY _raw_record_key
        """

    # Compatibility with snapshots created before compact raw-record IDs existed.
    return f"""
        SELECT
            _raw_record_key,
            row_number() OVER (ORDER BY _raw_record_key)::BIGINT AS _raw_record_id
        FROM (
            SELECT DISTINCT _raw_record_key
            FROM read_parquet({old_path})
        )
    """


def _duplicate_stats(records_path: Path) -> dict[str, Any]:
    con = duckdb.connect()
    try:
        duplicate_key_count, duplicate_row_count = con.execute(
            """
            WITH counts AS (
                SELECT _raw_record_key, count(*) AS n
                FROM read_parquet(?)
                GROUP BY _raw_record_key
            )
            SELECT
                count(*) FILTER (WHERE n > 1) AS duplicate_key_count,
                coalesce(sum(n - 1) FILTER (WHERE n > 1), 0) AS duplicate_row_count
            FROM counts
            """,
            [str(records_path)],
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
