"""Load frozen resolved Parquets into main's physical PostgreSQL contract.

DuckDB prepares dictionaries and relational COPY queries. PostgreSQL never
resolves an entity or reads a raw source body. Main's downstream functions own
any internal transactions; durable metadata records completed phase boundaries.
"""

from __future__ import annotations

from contextlib import closing
from dataclasses import asdict, dataclass, is_dataclass
import json
import logging
import os
import re
from pathlib import Path
from tempfile import TemporaryDirectory, gettempdir
from time import perf_counter

import duckdb
import psycopg2
from psycopg2 import sql
from psycopg2.extras import Json

from .locations import is_remote, measure_http_transfer

from .projection import companion_ddl, prepare_aligned_release
from .bulk import BulkCopyResult
from omnipath_subsets.db import validate_schema
from omnipath_subsets.runner import PRODUCTS, acquire_schema_lock, run_products
from .releases import load_release as read_release, verify_release
from .relational.db.schema import ensure_schema, ensure_source_partitions, ensure_deferred_indexes
from .relational.db.indexes import create_secondary_indexes

DIMENSION_IDS = {
    "data_source": "source_id",
    "dataset": "dataset_id",
    "vocab_identifier_type": "identifier_type_id",
    "vocab_entity_type": "entity_type_id",
    "vocab_relation_predicate": "relation_predicate_id",
    "vocab_relation_category": "relation_category_id",
    "vocab_annotation_scope": "annotation_scope_id",
    "vocab_entity_role": "entity_role_id",
    "vocab_resolution_status": "resolution_status_id",
}


@dataclass(frozen=True)
class LoadResult:
    release: str
    schema: str
    manifest_sha256: str
    resources: dict[str, str]
    counts: dict[str, int]
    phase_seconds: dict[str, float]
    committed_products: tuple[str, ...]
    compatibility: dict
    layout: str = "main"


def _emit(observer, event, **fields):
    if observer is not None:
        observer(event, **fields)


def _json(value):
    if is_dataclass(value):
        return _json(asdict(value))
    if isinstance(value, dict):
        return {str(k): _json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json(v) for v in value]
    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    return str(value)


def _lock(conn, schema):
    acquire_schema_lock(conn, schema)


def required_locks(partitioned: int, shared: int, sources: int) -> int:
    """Locks that deferring the constraints takes in one transaction.

    Every partitioned table and partitioned index gets a partition per source and
    a default one; dropping the constraints and indexes locks all of them, and the
    shared relations besides.
    """
    return partitioned * (sources + 1) + shared


def _check_lock_capacity(conn, schema, sources):
    """Fail before the projection when PostgreSQL cannot hold the deferral's locks.

    The shared lock table holds max_locks_per_transaction locks per connection slot
    (max_connections + max_prepared_transactions). With the default 64 and 100
    connections, 46 resources fit and 171 do not ("out of shared memory" after
    the projection, in _defer_constraints).
    """
    with conn.cursor() as cur:
        cur.execute(
            """SELECT count(*) FILTER (WHERE c.relkind IN ('p','I')),
                      count(*) FILTER (WHERE c.relkind IN ('r','i'))
            FROM pg_class c JOIN pg_namespace n ON n.oid=c.relnamespace
            WHERE n.nspname=%s AND NOT c.relispartition""",
            [schema],
        )
        partitioned, shared = cur.fetchone()
        cur.execute(
            "SELECT current_setting('max_locks_per_transaction')::int,"
            " current_setting('max_connections')::int"
            " + current_setting('max_prepared_transactions')::int"
        )
        per_slot, slots = cur.fetchone()
    needed = required_locks(partitioned, shared, sources)
    if needed > per_slot * slots:
        # A quarter for the foreign keys' referenced tables and other sessions.
        suggestion = 1 << max(6, (needed * 5 // 4 // slots).bit_length())
        raise ValueError(
            f"PostgreSQL's lock table holds {per_slot * slots} locks "
            f"(max_locks_per_transaction={per_slot} x {slots} slots); loading "
            f"{sources} resources needs about {needed}. Set max_locks_per_transaction "
            f"to {suggestion} or more and restart PostgreSQL."
        )


def _identity(conn, schema, expected=None):
    with conn.cursor() as cur:
        cur.execute(
            sql.SQL(
                "SELECT version,manifest_sha256 FROM {}.parquet_release WHERE singleton"
            ).format(sql.Identifier(schema))
        )
        row = cur.fetchone()
    if row is None or (expected is not None and tuple(row) != tuple(expected)):
        raise ValueError("Pinned PostgreSQL release identity changed")
    return tuple(row)


def _metadata_schema(conn, schema, release):
    with conn.cursor() as cur:
        cur.execute(
            sql.SQL("""CREATE TABLE {}.parquet_release (
            singleton boolean PRIMARY KEY DEFAULT true CHECK(singleton),
            version text NOT NULL, manifest_sha256 text NOT NULL,
            manifest jsonb NOT NULL, status text NOT NULL,
            counts jsonb NOT NULL DEFAULT '{{}}',
            phase_seconds jsonb NOT NULL DEFAULT '{{}}',
            compatibility jsonb NOT NULL DEFAULT '{{}}',
            constraint_plan jsonb NOT NULL DEFAULT '{{}}')""").format(sql.Identifier(schema))
        )
        cur.execute(
            sql.SQL("""CREATE TABLE {}.parquet_resource (
            resource text PRIMARY KEY, version text NOT NULL,
            manifest_sha256 text NOT NULL, manifest jsonb NOT NULL)""").format(
                sql.Identifier(schema)
            )
        )
        cur.execute(
            sql.SQL("""CREATE TABLE {}.parquet_phase (
            phase text PRIMARY KEY, completed_at timestamptz NOT NULL DEFAULT now(),
            seconds double precision NOT NULL, result jsonb NOT NULL)""").format(
                sql.Identifier(schema)
            )
        )
        cur.execute(
            sql.SQL(
                "INSERT INTO {}.parquet_release(singleton,version,manifest_sha256,manifest,status) VALUES(true,%s,%s,%s,%s)"
            ).format(sql.Identifier(schema)),
            [
                release.version,
                release.sha256,
                Json(json.loads(release.canonical_json)),
                "preparing",
            ],
        )
        cur.executemany(
            sql.SQL("INSERT INTO {}.parquet_resource VALUES(%s,%s,%s,%s)").format(
                sql.Identifier(schema)
            ),
            [
                (r.source, r.version, r.manifest_sha256, Json(json.loads(r.manifest_json)))
                for r in release.resources
            ],
        )
    conn.commit()


def _read_dimensions(conn, schema):
    result = {}
    with conn.cursor() as cur:
        for table, idcolumn in DIMENSION_IDS.items():
            columns = (idcolumn, "source_id", "name") if table == "dataset" else (idcolumn, "name")
            cur.execute(
                sql.SQL("SELECT {} FROM {}.{} ORDER BY {}").format(
                    sql.SQL(",").join(map(sql.Identifier, columns)),
                    sql.Identifier(schema),
                    sql.Identifier(table),
                    sql.Identifier(idcolumn),
                )
            )
            result[table] = tuple(cur.fetchall())
    return result


def _write_dimensions(conn, schema, dimensions):
    with conn.cursor() as cur:
        for table, rows in dimensions.items():
            idcolumn = DIMENSION_IDS[table]
            columns = (idcolumn, "source_id", "name") if table == "dataset" else (idcolumn, "name")
            if rows:
                cur.executemany(
                    sql.SQL("INSERT INTO {}.{}({}) VALUES({}) ON CONFLICT ({}) DO NOTHING").format(
                        sql.Identifier(schema),
                        sql.Identifier(table),
                        sql.SQL(",").join(map(sql.Identifier, columns)),
                        sql.SQL(",").join(sql.Placeholder() for _ in columns),
                        sql.Identifier(idcolumn),
                    ),
                    rows,
                )
            cur.execute("SELECT pg_get_serial_sequence(%s,%s)", [schema + "." + table, idcolumn])
            sequence = cur.fetchone()[0]
            if sequence and rows:
                cur.execute("SELECT setval(%s,%s,true)", [sequence, max(row[0] for row in rows)])
    conn.commit()


def _defer_constraints(conn, schema, tables):
    """Capture real main DDL, so deferred COPY restores the exact constraint contract."""
    with conn.cursor() as cur:
        cur.execute(
            """SELECT t.relname,c.conname,c.contype,pg_get_constraintdef(c.oid,true)
            FROM pg_constraint c JOIN pg_class t ON t.oid=c.conrelid
            JOIN pg_namespace n ON n.oid=t.relnamespace
            WHERE n.nspname=%s AND (t.relname=ANY(%s) OR c.contype='f') AND c.conparentid=0
                AND c.contype IN ('p','u','f','c','x')
            ORDER BY c.contype='f' DESC,t.relname,c.conname""",
            [schema, list(tables)],
        )
        constraints = cur.fetchall()
        for table, name, kind, definition in constraints:
            cur.execute(
                sql.SQL("ALTER TABLE {}.{} DROP CONSTRAINT {} CASCADE").format(
                    sql.Identifier(schema), sql.Identifier(table), sql.Identifier(name)
                )
            )
        cur.execute(
            """SELECT t.relname,i.relname,pg_get_indexdef(i.oid),t.relkind
            FROM pg_index x JOIN pg_class i ON i.oid=x.indexrelid
            JOIN pg_class t ON t.oid=x.indrelid JOIN pg_namespace n ON n.oid=t.relnamespace
            WHERE n.nspname=%s AND t.relname=ANY(%s)
              AND NOT EXISTS(SELECT 1 FROM pg_constraint c WHERE c.conindid=i.oid)
              AND NOT EXISTS(SELECT 1 FROM pg_inherits inh WHERE inh.inhrelid=i.oid)
            ORDER BY t.relname,i.relname""",
            [schema, list(tables)],
        )
        indexes = cur.fetchall()
        for table, name, definition, table_kind in indexes:
            cur.execute(
                sql.SQL("DROP INDEX {}.{} CASCADE").format(
                    sql.Identifier(schema), sql.Identifier(name)
                )
            )
        plan = {
            "constraints": constraints,
            "indexes": [(table, name, definition) for table, name, definition, _kind in indexes],
            "partitioned_index_definitions": [
                definition for _table, _name, definition, kind in indexes if kind == "p"
            ],
        }
        cur.execute(
            sql.SQL("UPDATE {}.parquet_release SET constraint_plan=%s WHERE singleton").format(
                sql.Identifier(schema)
            ),
            [Json(plan)],
        )
    conn.commit()
    return plan


def _restore_constraints(conn, schema, plan):
    with conn.cursor() as cur:
        for table, name, kind, definition in sorted(plan["constraints"], key=lambda r: r[2] == "f"):
            cur.execute(
                sql.SQL("ALTER TABLE {}.{} ADD CONSTRAINT {} ").format(
                    sql.Identifier(schema), sql.Identifier(table), sql.Identifier(name)
                )
                + sql.SQL(definition)
            )
        for table, name, definition in plan["indexes"]:
            if definition in plan.get("partitioned_index_definitions", ()):
                # pg_get_indexdef serializes a partitioned root as ON ONLY.
                # Main's original DDL builds a valid root and its child indexes.
                definition = definition.replace(" ON ONLY ", " ON ", 1)
            cur.execute(definition)
    # Validation and exact indexes are present before the base is advertised.


def _copy(conn, schema, item, con, spool):
    def literal(v):
        return "'" + str(v).replace("'", "''") + "'"

    with TemporaryDirectory(prefix=item.table + "-", dir=spool) as directory:
        output = Path(directory) / "csv"
        start = perf_counter()
        row = con.execute(
            f"COPY ({item.query}) TO {literal(output)} (FORMAT CSV, HEADER false, NULL '\\N', PER_THREAD_OUTPUT true, FILE_SIZE_BYTES '32MB', FILENAME_PATTERN 'part_{{i}}')"
        ).fetchone()
        if not row or type(row[0]) is not int or row[0] < 0:
            raise ValueError(f"DuckDB returned no valid staged row count for {item.table}")
        staged = row[0]
        stage_seconds = perf_counter() - start
        files = sorted(output.glob("*.csv"))
        if staged and not files:
            raise ValueError(f"DuckDB produced no CSV chunks for {item.table}")
        statement = sql.SQL("COPY {}.{} ({}) FROM STDIN WITH (FORMAT CSV, NULL '\\N')").format(
            sql.Identifier(schema),
            sql.Identifier(item.table),
            sql.SQL(",").join(map(sql.Identifier, item.columns)),
        )
        copied = size = 0
        start = perf_counter()
        for path in files:
            size += path.stat().st_size
            with conn.cursor() as cur, path.open("rb") as stream:
                cur.copy_expert(statement.as_string(conn), stream, size=1024 * 1024)
                if cur.rowcount < 0:
                    raise ValueError(f"PostgreSQL returned no COPY count for {item.table}")
                copied += cur.rowcount
        if copied != staged:
            raise ValueError(f"COPY count differs for {item.table}: {copied}/{staged}")
        return BulkCopyResult(staged, stage_seconds, perf_counter() - start, len(files), size)


def _ensure_molecular_type_index(conn, schema):
    """Restore main's exact molecular-type index, including an invalid root.

    Main's deferred helper omits this index. Older captured ON ONLY replay can
    leave its parent invalid with no attached children. Complete that state
    additively, preserving rows, the parent OID and any already attached work.
    """
    validate_schema(schema)
    index_name = "entity_evidence_resolution_molecular_type_idx"
    with conn.cursor() as cur:
        cur.execute(
            sql.SQL(
                "CREATE INDEX IF NOT EXISTS "
                "entity_evidence_resolution_molecular_type_idx "
                "ON {}.entity_evidence_resolution (molecular_type_id)"
            ).format(sql.Identifier(schema))
        )
        cur.execute(
            """SELECT ix.oid,i.indrelid,i.indisvalid,ix.relkind
            FROM pg_index i JOIN pg_class ix ON ix.oid=i.indexrelid
            JOIN pg_namespace ns ON ns.oid=ix.relnamespace
            WHERE ns.nspname=%s AND ix.relname=%s""",
            [schema, index_name],
        )
        parent_oid, table_oid, valid, index_kind = cur.fetchone()
        if not valid:
            if index_kind != "I":
                raise ValueError("Molecular-type index is invalid and not a partitioned root")
            cur.execute(
                """SELECT child.oid,ns.nspname,child.relname
                FROM pg_inherits inh JOIN pg_class child ON child.oid=inh.inhrelid
                JOIN pg_namespace ns ON ns.oid=child.relnamespace
                WHERE inh.inhparent=%s ORDER BY child.relname""",
                [table_oid],
            )
            partitions = cur.fetchall()
            for child_oid, child_schema, child_name in partitions:
                if child_schema != schema:
                    raise ValueError("Molecular-type source partition is outside the target schema")
                cur.execute(
                    """SELECT 1 FROM pg_inherits inh
                    JOIN pg_index attached ON attached.indexrelid=inh.inhrelid
                    WHERE inh.inhparent=%s AND attached.indrelid=%s""",
                    [parent_oid, child_oid],
                )
                if cur.fetchone() is not None:
                    continue
                matching = """SELECT ix.relname
                    FROM pg_index i JOIN pg_class ix ON ix.oid=i.indexrelid
                    JOIN pg_attribute col ON col.attrelid=i.indrelid
                         AND col.attname='molecular_type_id'
                    JOIN pg_index parent ON parent.indexrelid=%s
                    WHERE i.indrelid=%s AND i.indisvalid AND i.indisready
                      AND NOT i.indisunique AND i.indnkeyatts=1 AND i.indnatts=1
                      AND i.indkey[0]=col.attnum AND i.indclass=parent.indclass
                      AND i.indcollation=parent.indcollation AND i.indoption=parent.indoption
                      AND i.indpred IS NULL AND i.indexprs IS NULL
                      AND NOT EXISTS (SELECT 1 FROM pg_inherits inh WHERE inh.inhrelid=ix.oid)
                    ORDER BY ix.relname"""
                cur.execute(matching, [parent_oid, child_oid])
                candidates = cur.fetchall()
                if not candidates:
                    # Let PostgreSQL choose exactly the child name it would
                    # choose when building main's original partitioned index.
                    cur.execute(
                        sql.SQL("CREATE INDEX ON {}.{} (molecular_type_id)").format(
                            sql.Identifier(schema),
                            sql.Identifier(child_name),
                        )
                    )
                    cur.execute(matching, [parent_oid, child_oid])
                    candidates = cur.fetchall()
                if len(candidates) != 1:
                    raise ValueError("Molecular-type partition has no unique matching index")
                cur.execute(
                    sql.SQL("ALTER INDEX {}.{} ATTACH PARTITION {}.{}").format(
                        sql.Identifier(schema),
                        sql.Identifier(index_name),
                        sql.Identifier(schema),
                        sql.Identifier(candidates[0][0]),
                    )
                )
            cur.execute("SELECT indisvalid FROM pg_index WHERE indexrelid=%s", [parent_oid])
            if not cur.fetchone()[0]:
                raise ValueError("Molecular-type partitioned index remains invalid")
    conn.commit()


def _analyze(conn, schema):
    with conn.cursor() as cur:
        cur.execute(
            "SELECT c.relname FROM pg_class c JOIN pg_namespace n ON n.oid=c.relnamespace WHERE n.nspname=%s AND c.relkind IN ('r','p','m') ORDER BY c.relname",
            [schema],
        )
        for (table,) in cur.fetchall():
            cur.execute(
                sql.SQL("ANALYZE {}.{}").format(sql.Identifier(schema), sql.Identifier(table))
            )
    conn.commit()


def _drop_unfinished_schema(conn, schema):
    """Best effort: a failed cleanup must not hide the original load error."""
    try:
        conn.rollback()
        with conn.cursor() as cur:
            cur.execute("SELECT to_regclass(%s)", [f'"{schema}".parquet_phase'])
            if cur.fetchone()[0] is not None:
                cur.execute(
                    sql.SQL("SELECT 1 FROM {}.parquet_phase WHERE phase='base'").format(
                        sql.Identifier(schema)
                    )
                )
                if cur.fetchone() is not None:
                    return  # A committed base is resumable with finish.
            cur.execute(sql.SQL("DROP SCHEMA IF EXISTS {} CASCADE").format(sql.Identifier(schema)))
        conn.commit()
    except Exception:
        logging.getLogger(__name__).exception(
            "Could not drop unfinished schema %s; drop it manually before retrying", schema
        )


def _checkpoint(conn, schema, phase, seconds, result, identity, observer):
    _identity(conn, schema, identity)
    with conn.cursor() as cur:
        cur.execute(
            sql.SQL("""INSERT INTO {}.parquet_phase(phase,seconds,result) VALUES(%s,%s,%s)
            ON CONFLICT(phase) DO UPDATE SET seconds=EXCLUDED.seconds,result=EXCLUDED.result,completed_at=now()""").format(
                sql.Identifier(schema)
            ),
            [phase, seconds, Json(_json(result))],
        )
        cur.execute(
            sql.SQL("""UPDATE {}.parquet_release SET status=%s,
            phase_seconds=phase_seconds || jsonb_build_object(%s::text,%s::double precision) WHERE singleton""").format(
                sql.Identifier(schema)
            ),
            [phase, phase, seconds],
        )
    conn.commit()
    _emit(observer, "phase_committed", phase=phase, seconds=seconds, result=_json(result))


def finish_release(database_url, *, schema="omnipath", products=PRODUCTS, observer=None):
    """Resume committed main substrate; never repeat COPY or completed products."""
    validate_schema(schema)
    if not isinstance(products, (tuple, list)) or any(not isinstance(p, str) for p in products):
        raise ValueError("Products must be unique known products")
    if len(set(products)) != len(products) or any(p not in PRODUCTS for p in products):
        raise ValueError("Products must be unique known products")
    from .relational.pipeline import run_main_derivations

    with closing(psycopg2.connect(database_url)) as conn:
        _lock(conn, schema)
        identity = _identity(conn, schema)
        with conn.cursor() as cur:
            cur.execute(
                sql.SQL("SELECT phase FROM {}.parquet_phase").format(sql.Identifier(schema))
            )
            completed = {r[0] for r in cur.fetchall()}
            if "base" not in completed:
                raise ValueError("No committed base checkpoint; incomplete COPY cannot be resumed")
            cur.execute(
                sql.SQL(
                    "SELECT resource,manifest FROM {}.parquet_resource ORDER BY resource"
                ).format(sql.Identifier(schema))
            )
            metadata = dict(cur.fetchall())
        if "derived" not in completed:
            _emit(observer, "phase_start", phase="derived")
            start = perf_counter()
            outcome = run_main_derivations(
                conn, schema=schema, progress=True, published_metadata=metadata
            )
            _checkpoint(
                conn, schema, "derived", perf_counter() - start, outcome, identity, observer
            )
        run_products(
            database_url,
            schema,
            owner=conn,
            identity=identity,
            products=products,
            resume=True,
            observer=observer,
            progress=True,
        )
        with conn.cursor() as cur:
            cur.execute(
                sql.SQL("SELECT phase FROM {}.parquet_phase").format(sql.Identifier(schema))
            )
            completed = {r[0] for r in cur.fetchall()}
            if all(p in completed for p in PRODUCTS):
                cur.execute(
                    sql.SQL(
                        "UPDATE {}.parquet_release SET status='complete' WHERE singleton"
                    ).format(sql.Identifier(schema))
                )
            cur.execute(
                sql.SQL(
                    "SELECT version,manifest_sha256,counts,phase_seconds,compatibility FROM {}.parquet_release WHERE singleton"
                ).format(sql.Identifier(schema))
            )
            version, digest, counts, times, compat = cur.fetchone()
            cur.execute(
                sql.SQL("SELECT resource,version FROM {}.parquet_resource").format(
                    sql.Identifier(schema)
                )
            )
            resources = dict(cur.fetchall())
        conn.commit()
        return LoadResult(
            version,
            schema,
            digest,
            resources,
            counts,
            times,
            tuple(p for p in PRODUCTS if p in completed),
            compat,
        )


def load_release(
    data_root,
    manifest_path,
    database_url,
    *,
    schema="omnipath",
    duckdb_threads=1,
    memory_limit="1GB",
    temp_directory=None,
    defer_constraints=True,
    base_only=False,
    products=PRODUCTS,
    observer=None,
    retain_published_provenance=False,
):
    """Create main tables and required molecular/quantity context, then finish a validated base.

    Exact published audit copies are opt-in; all original details remain in the
    immutable release Parquets. Existing checkpoints keep their stored layout
    and counts when passed to :func:`finish_release`.
    """
    validate_schema(schema)
    if type(duckdb_threads) is not int or duckdb_threads < 1:
        raise ValueError("duckdb_threads must be a positive integer")
    if not isinstance(memory_limit, str) or not re.fullmatch(
        r"[1-9][0-9]*(?:\.[0-9]+)?\s*(?:KB|MB|GB|KiB|MiB|GiB)", memory_limit
    ):
        raise ValueError("memory_limit must be a positive size such as 512MB")
    if type(defer_constraints) is not bool or type(base_only) is not bool:
        raise ValueError("defer_constraints and base_only must be booleans")
    if type(retain_published_provenance) is not bool:
        raise ValueError("retain_published_provenance must be a boolean")
    if not isinstance(products, (tuple, list)) or any(not isinstance(p, str) for p in products):
        raise ValueError("Products must be unique known products")
    if len(set(products)) != len(products) or any(p not in PRODUCTS for p in products):
        raise ValueError("Products must be unique known products")
    start = perf_counter()
    with measure_http_transfer() as initial_transfer:
        release = read_release(data_root, manifest_path)
    artifact_seconds = perf_counter() - start
    _emit(observer, "artifact_validation_complete", **asdict(initial_transfer))
    staging_parent = (
        Path(manifest_path).parent
        if not is_remote(manifest_path) and os.access(Path(manifest_path).parent, os.W_OK)
        else Path(gettempdir())
    )
    spool = Path(temp_directory or staging_parent / "main-staging")
    spool.mkdir(parents=True, exist_ok=True)
    with closing(psycopg2.connect(database_url)) as conn:
        _lock(conn, schema)
        with conn.cursor() as cur:
            cur.execute("SELECT EXISTS(SELECT 1 FROM pg_namespace WHERE nspname=%s)", [schema])
            if cur.fetchone()[0]:
                raise ValueError(
                    "Aligned loading requires a new schema; use finish for a checkpoint"
                )
        try:
            ensure_schema(conn, schema=schema, indexes=False, progress=True)
            _check_lock_capacity(conn, schema, len(release.resources))
            _metadata_schema(conn, schema, release)
            identity = (release.version, release.sha256)
            with TemporaryDirectory(prefix="main-projection-", dir=spool) as directory:
                with duckdb.connect(str(Path(directory) / "projection.duckdb")) as con:
                    con.execute("SET threads=?", [duckdb_threads])
                    con.execute("SET memory_limit=?", [memory_limit])
                    con.execute("SET temp_directory=?", [str(Path(directory) / "spill")])
                    _emit(observer, "phase_start", phase="projection")
                    start = perf_counter()
                    plan = prepare_aligned_release(
                        con,
                        release,
                        dimension_rows=_read_dimensions(conn, schema),
                        progress=lambda fields: _emit(observer, "projection_progress", **fields),
                        retain_published_provenance=retain_published_provenance,
                    )
                    projection_seconds = perf_counter() - start
                    _write_dimensions(conn, schema, plan.dimensions)
                    for _, source in plan.sources:
                        ensure_source_partitions(conn, schema=schema, source=source)
                    with conn.cursor() as cur:
                        for statement in companion_ddl(
                            schema,
                            retain_published_provenance=retain_published_provenance,
                        ):
                            cur.execute(statement)
                    conn.commit()
                    constraint_plan = (
                        _defer_constraints(conn, schema, [q.table for q in plan.queries])
                        if defer_constraints
                        else {"constraints": [], "indexes": []}
                    )
                    start = perf_counter()
                    counts = {}
                    staging = copy_seconds = 0.0
                    for item in plan.queries:
                        _emit(
                            observer, "copy_start", table=item.table, rows=plan.counts[item.table]
                        )
                        copied = _copy(conn, schema, item, con, spool)
                        if copied.rows != plan.counts[item.table]:
                            raise ValueError("Prepared projection changed during COPY")
                        counts[item.table] = copied.rows
                        staging += copied.stage_seconds
                        copy_seconds += copied.copy_seconds
                        _emit(observer, "copy_complete", table=item.table, **asdict(copied))
                    _emit(observer, "phase_start", phase="constraints")
                    constraints_start = perf_counter()
                    _restore_constraints(conn, schema, constraint_plan)
                    constraints_seconds = perf_counter() - constraints_start
                    # Ensure helpers can commit only after every captured key/FK was restored.
                    conn.commit()
                    ensure_deferred_indexes(conn, schema=schema, progress=True)
                    _ensure_molecular_type_index(conn, schema)
                    create_secondary_indexes(conn, schema=schema)
                    _analyze(conn, schema)
                    with measure_http_transfer() as final_transfer:
                        verify_release(release)
                    _emit(observer, "artifact_revalidation_complete", **asdict(final_transfer))
                    with conn.cursor() as cur:
                        cur.execute(
                            sql.SQL(
                                "UPDATE {}.parquet_release SET counts=%s,compatibility=%s,phase_seconds=%s WHERE singleton"
                            ).format(sql.Identifier(schema)),
                            [
                                Json(counts),
                                Json(dict(plan.compatibility)),
                                Json(
                                    {
                                        "artifact_validation": artifact_seconds,
                                        "artifact_revalidation": final_transfer.seconds,
                                        "projection": projection_seconds,
                                        "csv_staging": staging,
                                        "copy": copy_seconds,
                                        "constraints": constraints_seconds,
                                    }
                                ),
                            ],
                        )
                    _checkpoint(
                        conn,
                        schema,
                        "base",
                        perf_counter() - start,
                        {"counts": counts},
                        identity,
                        observer,
                    )
        except BaseException:
            # Nothing durable exists before the base checkpoint, and this call
            # created the schema itself under the lock: remove it so a rerun works.
            _drop_unfinished_schema(conn, schema)
            raise

    if base_only:
        return LoadResult(
            release.version,
            schema,
            release.sha256,
            {r.source: r.version for r in release.resources},
            counts,
            {
                "artifact_validation": artifact_seconds,
                "artifact_revalidation": final_transfer.seconds,
                "projection": projection_seconds,
                "csv_staging": staging,
                "copy": copy_seconds,
                "constraints": constraints_seconds,
            },
            (),
            dict(plan.compatibility),
        )
    return finish_release(database_url, schema=schema, products=products, observer=observer)
