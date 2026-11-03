"""Aligned loader rejects invalid requests and preserves COPY failure boundaries."""

import csv
import io
from contextlib import closing
import uuid
import json

import duckdb
import pytest
import psycopg2
from psycopg2 import sql

from omnipath_postgres import loader as aligned_loader
from omnipath_subsets import runner
from omnipath_postgres.projection import CopyQuery


def forbid_io(*args, **kwargs):
    pytest.fail("Invalid loader parameters must not read artifacts or connect")


@pytest.mark.parametrize("schema", ["public", "pg_temp", "", "bad-name", None])
def test_invalid_schema_before_io(monkeypatch, schema):
    monkeypatch.setattr(aligned_loader, "read_release", forbid_io)
    monkeypatch.setattr(aligned_loader.psycopg2, "connect", forbid_io)
    with pytest.raises(ValueError, match="schema"):
        aligned_loader.load_release("missing", "missing", "unused", schema=schema)


@pytest.mark.parametrize("threads", [0, -1, True, 1.5, None])
def test_invalid_duckdb_threads_before_io(monkeypatch, threads):
    monkeypatch.setattr(aligned_loader, "read_release", forbid_io)
    monkeypatch.setattr(aligned_loader.psycopg2, "connect", forbid_io)
    with pytest.raises(ValueError, match="duckdb_threads"):
        aligned_loader.load_release("missing", "missing", "unused", duckdb_threads=threads)


@pytest.mark.parametrize(
    "limit", [None, "", 0, True, "0MB", "-1MB", "unlimited", "512MB; SELECT 1"]
)
def test_invalid_memory_limit_before_io(monkeypatch, limit):
    monkeypatch.setattr(aligned_loader, "read_release", forbid_io)
    monkeypatch.setattr(aligned_loader.psycopg2, "connect", forbid_io)
    with pytest.raises(ValueError, match="memory_limit"):
        aligned_loader.load_release("missing", "missing", "unused", memory_limit=limit)


@pytest.mark.parametrize("field", ["defer_constraints", "base_only", "retain_published_provenance"])
@pytest.mark.parametrize("value", [None, 0, 1, "true"])
def test_nonboolean_modes_before_io(monkeypatch, field, value):
    monkeypatch.setattr(aligned_loader, "read_release", forbid_io)
    monkeypatch.setattr(aligned_loader.psycopg2, "connect", forbid_io)
    with pytest.raises(ValueError, match=field):
        aligned_loader.load_release("missing", "missing", "unused", **{field: value})


@pytest.mark.parametrize("products", [("cosmos", "cosmos"), ("unknown",), "cosmos", None])
def test_invalid_products_before_io(monkeypatch, products):
    monkeypatch.setattr(aligned_loader, "read_release", forbid_io)
    monkeypatch.setattr(aligned_loader.psycopg2, "connect", forbid_io)
    with pytest.raises(ValueError, match="[Pp]roducts"):
        aligned_loader.load_release("missing", "missing", "unused", products=products)


@pytest.mark.parametrize("products", [("cosmos", "cosmos"), ("unknown",), "cosmos", None])
def test_invalid_finish_products_before_connect(monkeypatch, products):
    monkeypatch.setattr(aligned_loader.psycopg2, "connect", forbid_io)
    with pytest.raises(ValueError, match="[Pp]roducts"):
        aligned_loader.finish_release("unused", schema="fixture", products=products)


class Cursor:
    def __init__(self, connection):
        self.connection = connection
        self.rowcount = -1

    def __enter__(self):
        return self

    def __exit__(self, *_):
        pass

    def copy_expert(self, statement, stream, **kwargs):
        self.connection.statements.append(statement)
        payload = stream.read()
        self.connection.payloads.append(payload)
        if self.connection.failure:
            raise RuntimeError("deliberate PostgreSQL COPY failure")
        parsed = list(csv.reader(io.StringIO(payload.decode())))
        self.rowcount = (
            self.connection.rowcounts.pop(0) if self.connection.rowcounts else len(parsed)
        )


class Connection:
    def __init__(self, *, rowcounts=(), failure=False):
        self.rowcounts = list(rowcounts)
        self.failure = failure
        self.payloads = []
        self.statements = []
        self.commits = 0

    def cursor(self):
        return Cursor(self)

    def commit(self):
        self.commits += 1


@pytest.fixture
def fake_copy_render(monkeypatch):
    # Only transport is mocked. DuckDB stages real CSV. Real PostgreSQL DDL,
    # identifiers and driver COPY are independently covered by the server pilot.
    monkeypatch.setattr(aligned_loader.sql.Composed, "as_string", lambda self, conn: "COPY fixture")


def test_zero_rows_copy_is_empty_success_and_no_hidden_commit(tmp_path, fake_copy_render):
    conn = Connection()
    with duckdb.connect() as con:
        item = CopyQuery("entity", ("entity_id",), "SELECT NULL::UUID WHERE false")
        result = aligned_loader._copy(conn, "fixture", item, con, tmp_path)
    assert result.rows == 0
    assert conn.commits == 0
    assert list(tmp_path.iterdir()) == []


def test_copy_null_literal_null_token_quotes_and_unicode_are_distinct(tmp_path, fake_copy_render):
    conn = Connection()
    with duckdb.connect() as con:
        con.execute("CREATE TABLE fixture(value VARCHAR)")
        con.executemany(
            "INSERT INTO fixture VALUES (?)", [(None,), (r"\N",), ("α,'\"\nnewline",), ("",)]
        )
        item = CopyQuery("annotation", ("value",), "SELECT value FROM fixture")
        result = aligned_loader._copy(conn, "fixture", item, con, tmp_path)
    assert result.rows == 4
    payload = b"".join(conn.payloads).decode()
    assert "\\N\n" in payload
    assert '"\\N"\n' in payload
    assert '""\n' in payload
    assert "α" in payload
    assert conn.commits == 0
    assert list(tmp_path.iterdir()) == []


def test_copy_count_mismatch_cleans_staging_and_does_not_commit(tmp_path, fake_copy_render):
    conn = Connection(rowcounts=(2,))
    with duckdb.connect() as con:
        item = CopyQuery("annotation", ("value",), "SELECT 'one'")
        with pytest.raises(ValueError, match="COPY count differs"):
            aligned_loader._copy(conn, "fixture", item, con, tmp_path)
    assert conn.commits == 0
    assert list(tmp_path.iterdir()) == []


def test_copy_transport_failure_cleans_staging_and_does_not_commit(tmp_path, fake_copy_render):
    conn = Connection(failure=True)
    with duckdb.connect() as con:
        item = CopyQuery("annotation", ("value",), "SELECT 'one'")
        with pytest.raises(RuntimeError, match="deliberate PostgreSQL COPY failure"):
            aligned_loader._copy(conn, "fixture", item, con, tmp_path)
    assert conn.commits == 0
    assert list(tmp_path.iterdir()) == []


class InvalidStaging:
    def __init__(self, result):
        self.result = result

    def execute(self, _):
        return self

    def fetchone(self):
        return self.result


@pytest.mark.parametrize("result", [None, (), (-1,), (True,), ("1",)])
def test_invalid_duckdb_stage_count_is_rejected(tmp_path, fake_copy_render, result):
    conn = Connection()
    item = CopyQuery("annotation", ("value",), "SELECT 'one'")
    with pytest.raises(ValueError, match="staged row count"):
        aligned_loader._copy(conn, "fixture", item, InvalidStaging(result), tmp_path)
    assert not conn.statements
    assert conn.commits == 0
    assert list(tmp_path.iterdir()) == []


def test_nonempty_staged_count_requires_csv_chunks(tmp_path, fake_copy_render):
    conn = Connection()
    item = CopyQuery("annotation", ("value",), "SELECT 'one'")
    with pytest.raises(ValueError, match="CSV chunks"):
        aligned_loader._copy(conn, "fixture", item, InvalidStaging((1,)), tmp_path)
    assert list(tmp_path.iterdir()) == []


@pytest.fixture
def checkpoint_database(postgres_dsn):
    """Tiny disposable publication schema; no resource build or source data."""
    schema = "checkpoint_" + uuid.uuid4().hex
    with closing(psycopg2.connect(postgres_dsn)) as connection:
        with connection.cursor() as cursor:
            cursor.execute(sql.SQL("CREATE SCHEMA {}").format(sql.Identifier(schema)))
            cursor.execute(
                sql.SQL("""CREATE TABLE {}.parquet_release (
                singleton boolean PRIMARY KEY, version text,manifest_sha256 text,
                status text,phase_seconds jsonb NOT NULL DEFAULT '{{}}')""").format(
                    sql.Identifier(schema)
                )
            )
            cursor.execute(
                sql.SQL("""CREATE TABLE {}.parquet_phase (
                phase text PRIMARY KEY,seconds double precision,result jsonb,
                completed_at timestamptz DEFAULT now())""").format(sql.Identifier(schema))
            )
            cursor.execute(
                sql.SQL("CREATE TABLE {}.product_contents(value text)").format(
                    sql.Identifier(schema)
                )
            )
            cursor.execute(
                sql.SQL(
                    "INSERT INTO {}.parquet_release VALUES(true,'release','digest','derived','{{}}')"
                ).format(sql.Identifier(schema))
            )
            cursor.execute(
                sql.SQL("INSERT INTO {}.parquet_phase(phase) VALUES('base'),('derived')").format(
                    sql.Identifier(schema)
                )
            )
        connection.commit()
        try:
            yield postgres_dsn, schema, connection
        finally:
            connection.rollback()
            with connection.cursor() as cursor:
                cursor.execute(sql.SQL("DROP SCHEMA {} CASCADE").format(sql.Identifier(schema)))
            connection.commit()


@pytest.mark.parametrize("wrapped", [True, False])
@pytest.mark.parametrize("observer_fails", [False, True])
def test_product_checkpoint_commits_tables_and_both_metadata_interfaces_before_observer(
    checkpoint_database,
    wrapped,
    observer_fails,
    monkeypatch,
):
    dsn, schema, connection = checkpoint_database
    stats = {"rows": 1, "sources": ["fixture"]}
    outcome = {"phase_seconds": {"projection": 0.25}, "result": stats} if wrapped else stats
    outcome = outcome if wrapped else {"result": stats}

    def build(worker, *args, **kwargs):
        with worker.cursor() as cursor:
            cursor.execute(
                sql.SQL("INSERT INTO {}.product_contents VALUES('complete product')").format(
                    sql.Identifier(schema)
                )
            )
        return outcome

    monkeypatch.setattr(runner, "run_product", build)
    times = iter((0.0, 1.5))
    monkeypatch.setattr(runner, "perf_counter", lambda: next(times))
    runner.acquire_schema_lock(connection, schema)
    observed = []

    def observer(event, **fields):
        if event == "phase_start":
            return
        # A separate connection can see only actually committed publication.
        with closing(psycopg2.connect(dsn)) as reader, reader.cursor() as cursor:
            cursor.execute(
                sql.SQL("SELECT value FROM {}.product_contents").format(sql.Identifier(schema))
            )
            assert cursor.fetchall() == [("complete product",)]
            cursor.execute(
                sql.SQL(
                    "SELECT product,release_id,manifest_sha256,stats FROM {}.subset_build_metadata"
                ).format(sql.Identifier(schema))
            )
            assert cursor.fetchall() == [("cosmos", "release", "digest", stats)]
            cursor.execute(
                sql.SQL(
                    "SELECT phase,seconds,result FROM {}.parquet_phase WHERE phase='cosmos'"
                ).format(sql.Identifier(schema))
            )
            assert cursor.fetchall() == [("cosmos", 1.5, outcome)]
            cursor.execute(
                sql.SQL("SELECT status,phase_seconds FROM {}.parquet_release").format(
                    sql.Identifier(schema)
                )
            )
            assert cursor.fetchall() == [("cosmos", {"cosmos": 1.5})]
        observed.append((event, fields))
        if observer_fails:
            raise OSError("deliberate observer failure")

    if observer_fails:
        with pytest.raises(OSError, match="observer failure"):
            runner.run_products(
                dsn,
                schema,
                owner=connection,
                identity=("release", "digest"),
                products=("cosmos",),
                observer=observer,
            )
    else:
        runner.run_products(
            dsn,
            schema,
            owner=connection,
            identity=("release", "digest"),
            products=("cosmos",),
            observer=observer,
        )
    assert observed == [("phase_committed", {"phase": "cosmos", "seconds": 1.5, "result": outcome})]
    # Caller rollback after an observer exception must retain published work.
    connection.rollback()
    with connection.cursor() as cursor:
        cursor.execute(
            sql.SQL("SELECT count(*) FROM {}.product_contents").format(sql.Identifier(schema))
        )
        assert cursor.fetchone()[0] == 1


def test_checkpoint_identity_change_never_commits_pending_product(checkpoint_database):
    dsn, schema, connection = checkpoint_database
    with connection.cursor() as cursor:
        cursor.execute(
            sql.SQL("INSERT INTO {}.product_contents VALUES('unfinished product')").format(
                sql.Identifier(schema)
            )
        )
    with pytest.raises(ValueError, match="release changed"):
        runner._publish(
            connection, schema, "cosmos", ("wrong release", "wrong digest"), {"result": {}}, 1.5
        )
    with closing(psycopg2.connect(dsn)) as reader, reader.cursor() as cursor:
        cursor.execute(
            sql.SQL("SELECT count(*) FROM {}.product_contents").format(sql.Identifier(schema))
        )
        assert cursor.fetchone()[0] == 0
        cursor.execute(
            sql.SQL("SELECT count(*) FROM {}.parquet_phase WHERE phase='cosmos'").format(
                sql.Identifier(schema)
            )
        )
        assert cursor.fetchone()[0] == 0


def test_current_clis_finish_without_parquets_skip_committed_then_explicitly_rebuild(
    checkpoint_database,
    monkeypatch,
    capsys,
):
    """Both commands publish through the runner against a tiny pinned substrate."""
    from omnipath_postgres import cli as postgres_cli
    from omnipath_subsets import cli as subsets_cli

    dsn, schema, connection = checkpoint_database
    with connection.cursor() as cursor:
        cursor.execute(
            sql.SQL(
                "ALTER TABLE {}.parquet_release ADD COLUMN counts jsonb DEFAULT '{{}}', ADD COLUMN compatibility jsonb DEFAULT '{{}}'"
            ).format(sql.Identifier(schema))
        )
        cursor.execute(
            sql.SQL("UPDATE {}.parquet_release SET version='2026.09'").format(
                sql.Identifier(schema)
            )
        )
        cursor.execute(
            sql.SQL(
                "CREATE TABLE {}.parquet_resource(resource text,version text,manifest jsonb)"
            ).format(sql.Identifier(schema))
        )
        cursor.execute(
            sql.SQL("INSERT INTO {}.parquet_resource VALUES ('fixture','1.0.0','{{}}')").format(
                sql.Identifier(schema)
            )
        )
    connection.commit()
    monkeypatch.setattr(aligned_loader, "read_release", forbid_io)
    monkeypatch.setattr(aligned_loader, "_copy", forbid_io)

    def product(worker, product, **kwargs):
        assert product == "cosmos"
        with worker.cursor() as cursor:
            cursor.execute(
                sql.SQL("INSERT INTO {}.product_contents VALUES('complete product')").format(
                    sql.Identifier(schema)
                )
            )
            cursor.execute(
                sql.SQL("SELECT count(*) FROM {}.product_contents").format(sql.Identifier(schema))
            )
            rows = cursor.fetchone()[0]
        return {"result": {"rows": rows}, "phase_seconds": 0.0}

    monkeypatch.setattr(runner, "run_product", product)
    args = ["--database-url", dsn, "--schema", schema, "--products", "cosmos"]
    assert postgres_cli.main(["--finish", *args]) == 0
    first = json.loads(capsys.readouterr().out)
    assert first["release"] == "2026.09" and first["resources"] == {"fixture": "1.0.0"}
    assert first["committed_products"] == ["cosmos"]
    assert postgres_cli.main(["--finish", *args]) == 0
    second = json.loads(capsys.readouterr().out)
    assert second["phase_seconds"] == first["phase_seconds"]
    assert subsets_cli.main(["build", *args]) == 0
    rebuilt = json.loads(capsys.readouterr().out)
    assert rebuilt["release"] == "2026.09" and rebuilt["products"] == {"cosmos": {"rows": 2}}
    with connection.cursor() as cursor:
        cursor.execute(
            sql.SQL("SELECT count(*) FROM {}.product_contents").format(sql.Identifier(schema))
        )
        assert cursor.fetchone()[0] == 2
        cursor.execute(
            sql.SQL(
                "SELECT release_id,manifest_sha256,stats FROM {}.subset_build_metadata WHERE product='cosmos'"
            ).format(sql.Identifier(schema))
        )
        assert cursor.fetchone() == ("2026.09", "digest", {"rows": 2})


class _LockSettingsCursor:
    """Answers the two queries of _check_lock_capacity."""

    def __init__(self, relations, settings):
        self.answers = [relations, settings]

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def execute(self, *args):
        pass

    def fetchone(self):
        return self.answers.pop(0)


class _LockSettingsConnection:
    def __init__(self, relations, settings):
        self.relations, self.settings = relations, settings

    def cursor(self):
        return _LockSettingsCursor(self.relations, self.settings)


def test_lock_capacity_is_checked_before_the_projection():
    # dev5, October 2026: 59 partitioned tables and indexes, ~290 shared relations;
    # the default 64 locks x 100 connections held 46 resources but not 171.
    assert aligned_loader.required_locks(59, 290, 46) < 64 * 100 < aligned_loader.required_locks(59, 290, 171)
    aligned_loader._check_lock_capacity(_LockSettingsConnection((59, 290), (64, 100)), "omnipath", 46)
    with pytest.raises(ValueError, match="max_locks_per_transaction to 256"):
        aligned_loader._check_lock_capacity(_LockSettingsConnection((59, 290), (64, 100)), "omnipath", 171)
    aligned_loader._check_lock_capacity(_LockSettingsConnection((59, 290), (4096, 100)), "omnipath", 171)
