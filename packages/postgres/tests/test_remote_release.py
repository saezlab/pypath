"""Remote pins retain validation and full relational projection parity."""

from contextlib import contextmanager
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
import json
from pathlib import Path
from threading import Thread

import duckdb
import pytest

from omnipath_postgres.projection import prepare_aligned_release
from omnipath_postgres.locations import HTTPRangeFile, validate_url
from omnipath_postgres.releases import FILES, ReleaseValidationError, load_release, verify_release
from omnipath_core.fixtures import write_manifest
from published_fixture import write_fixture


@contextmanager
def endpoint(root, *, ranges=True):
    class Handler(SimpleHTTPRequestHandler):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, directory=str(root), **kwargs)

        def log_message(self, *args):
            pass

        def do_GET(self):
            value = self.headers.get("Range")
            if not value or not ranges:
                return super().do_GET()
            path = Path(self.translate_path(self.path))
            if not path.is_file():
                return self.send_error(404)
            start, end = map(int, value.removeprefix("bytes=").split("-"))
            size = path.stat().st_size
            end = min(end, size - 1)
            self.send_response(206)
            self.send_header("Content-Length", str(end - start + 1))
            self.send_header("Content-Range", f"bytes {start}-{end}/{size}")
            self.end_headers()
            with path.open("rb") as stream:
                stream.seek(start)
                self.wfile.write(stream.read(end - start + 1))

    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    worker = Thread(target=server.serve_forever, daemon=True)
    worker.start()
    try:
        yield f"http://127.0.0.1:{server.server_port}"
    finally:
        server.shutdown()
        server.server_close()
        worker.join()


@pytest.fixture
def published(tmp_path):
    directory = tmp_path / "resources/fixture/1"
    directory.mkdir(parents=True)
    write_fixture(directory)
    write_manifest(directory, "fixture", "1")
    release = tmp_path / "releases/2026.10.json"
    release.parent.mkdir()
    release.write_text(
        json.dumps({"schema_version": 1, "version": "2026.10", "resources": {"fixture": "1"}})
    )
    return tmp_path, release, directory


def test_remote_validation_and_projection_match_local_complete_rows(published):
    root, path, _ = published
    local = load_release(root, path)
    with endpoint(root) as base:
        remote = load_release(base, base + "/releases/2026.10.json")
        assert remote.sha256 == local.sha256
        assert remote.resources[0].manifest_sha256 == local.resources[0].manifest_sha256
        assert remote.resources[0].table_path("entity").startswith(base)
        verify_release(remote)
        with duckdb.connect() as left, duckdb.connect() as right:
            a = prepare_aligned_release(left, local)
            b = prepare_aligned_release(right, remote)
            assert a.counts == b.counts
            assert a.dimensions == b.dimensions
            assert a.compatibility == b.compatibility
            for x, y in zip(a.queries, b.queries, strict=True):
                assert x.table == y.table and x.columns == y.columns
                assert sorted(map(repr, left.execute(x.query).fetchall())) == sorted(
                    map(repr, right.execute(y.query).fetchall())
                )
    # No Parquet download/cache has been created alongside the immutable files.
    assert sorted(p.name for p in root.iterdir()) == ["releases", "resources"]


def test_local_manifest_can_pin_remote_artifacts(published):
    root, path, _ = published
    with endpoint(root) as base:
        remote = load_release(base, path)
        assert remote.manifest_path == path
        verify_release(remote)


def test_remote_mutation_is_rejected_before_checkpoint(published):
    root, path, directory = published
    with endpoint(root) as base:
        remote = load_release(base, path)
        target = directory / "evidence_payloads.parquet"
        data = bytearray(target.read_bytes())
        data[10] ^= 1
        target.write_bytes(data)
        with pytest.raises(ReleaseValidationError):
            verify_release(remote)


def test_remote_manifest_change_is_rejected(published):
    root, path, _ = published
    with endpoint(root) as base:
        remote = load_release(base, base + "/releases/2026.10.json")
        path.write_text(path.read_text() + "\n")
        with pytest.raises(ReleaseValidationError, match="manifest changed"):
            verify_release(remote)


def test_missing_range_support_is_rejected(published):
    root, _, _ = published
    with endpoint(root, ranges=False) as base:
        with HTTPRangeFile(base + "/resources/fixture/1/entity.parquet") as stream:
            with pytest.raises(OSError, match="byte ranges"):
                stream.read(4)


@pytest.mark.parametrize(
    "url",
    [
        "http://example.org/data",
        "https://user:secret@example.org/data",
        "https://example.org/data?latest=true",
        "https://example.org/data#latest",
        "https://example.org/data/../other",
        "https://example.org/data/%2e%2e/other",
    ],
)
def test_remote_urls_are_literal_and_credential_free(url):
    with pytest.raises(ValueError):
        validate_url(url)


@pytest.mark.integration
def test_remote_base_load_commits_actual_postgres_rows(published, postgres_dsn):
    from contextlib import closing
    import uuid
    import psycopg2
    from psycopg2 import sql
    from omnipath_postgres.loader import load_release

    root, _, _ = published
    schema = "remote_fixture_" + uuid.uuid4().hex
    with closing(psycopg2.connect(postgres_dsn)) as conn:
        with conn.cursor() as cur:
            cur.execute("CREATE EXTENSION IF NOT EXISTS pg_trgm")
        conn.commit()
        try:
            with endpoint(root) as base:
                result = load_release(
                    base,
                    base + "/releases/2026.10.json",
                    postgres_dsn,
                    schema=schema,
                    base_only=True,
                    temp_directory=root / "spool",
                )
            with conn.cursor() as cur:
                for table, expected in result.counts.items():
                    cur.execute(
                        sql.SQL("SELECT count(*) FROM {}.{}").format(
                            sql.Identifier(schema), sql.Identifier(table)
                        )
                    )
                    assert cur.fetchone()[0] == expected
                cur.execute(
                    sql.SQL("SELECT phase FROM {}.parquet_phase").format(sql.Identifier(schema))
                )
                assert cur.fetchall() == [("base",)]
                cur.execute(
                    sql.SQL("SELECT manifest FROM {}.parquet_release WHERE singleton").format(
                        sql.Identifier(schema)
                    )
                )
                assert cur.fetchone()[0]["resources"] == {"fixture": "1"}
        finally:
            conn.rollback()
            with conn.cursor() as cur:
                cur.execute(
                    sql.SQL("DROP SCHEMA IF EXISTS {} CASCADE").format(sql.Identifier(schema))
                )
            conn.commit()


def test_validation_transfer_metrics_are_actual_bytes_and_scope_bound(published):
    from omnipath_postgres.locations import measure_http_transfer, read_bytes

    root, release, directory = published
    total_artifact_bytes = sum((directory / name).stat().st_size for name in FILES)
    with endpoint(root) as base:
        with measure_http_transfer() as initial:
            pinned = load_release(base, base + "/releases/2026.10.json")
        with measure_http_transfer() as final:
            verify_release(pinned)
        assert initial.requests > len(FILES)
        assert initial.bytes_read >= total_artifact_bytes + release.stat().st_size
        assert final.bytes_read == initial.bytes_read
        assert initial.seconds > 0 and final.seconds > 0
        before = initial.bytes_read
        assert read_bytes(base + "/releases/2026.10.json") == release.read_bytes()
        assert initial.bytes_read == before


def test_bounded_measurement_cli_reports_projection_and_validation(published, capsys, monkeypatch):
    import importlib.util
    import sys

    root, release, _ = published
    script = Path(__file__).resolve().parents[3] / "scripts/measure_parquet_input.py"
    spec = importlib.util.spec_from_file_location("measure_input", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    with endpoint(root) as base:
        monkeypatch.setattr(
            sys, "argv", [str(script), base + "/releases/2026.10.json", "--data-root", base]
        )
        module.main()
    result = json.loads(capsys.readouterr().out)
    assert result["initial_validation"]["bytes_read"] >= result["artifact_bytes"]
    assert result["final_validation"]["bytes_read"] > 0
    assert result["projection_seconds"] > 0
    assert result["projected_counts"]["entity"] > 0
    monkeypatch.setattr(
        sys, "argv", [str(script), str(release), "--data-root", str(root), "--max-bytes", "1"]
    )
    with pytest.raises(SystemExit):
        module.main()


def test_failed_base_load_removes_its_schema_so_a_rerun_succeeds(
    published, postgres_dsn, monkeypatch
):
    from contextlib import closing
    import uuid
    import psycopg2
    from omnipath_postgres import loader

    root, _, _ = published
    schema = "failed_fixture_" + uuid.uuid4().hex

    def exists():
        with closing(psycopg2.connect(postgres_dsn)) as conn, conn.cursor() as cur:
            cur.execute("SELECT EXISTS(SELECT 1 FROM pg_namespace WHERE nspname=%s)", [schema])
            return cur.fetchone()[0]

    def fail(*args, **kwargs):
        raise RuntimeError("deliberate constraint failure")

    arguments = (root, root / "releases/2026.10.json", postgres_dsn)
    options = dict(schema=schema, base_only=True, temp_directory=root / "spool")
    with monkeypatch.context() as patch:
        patch.setattr(loader, "_restore_constraints", fail)
        with pytest.raises(RuntimeError, match="deliberate"):
            loader.load_release(*arguments, **options)
    assert not exists()
    try:
        assert loader.load_release(*arguments, **options).counts
    finally:
        with closing(psycopg2.connect(postgres_dsn)) as conn, conn.cursor() as cur:
            cur.execute(f'DROP SCHEMA IF EXISTS "{schema}" CASCADE')
            conn.commit()
