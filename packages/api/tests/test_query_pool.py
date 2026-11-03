"""Resource limits apply across request threads and survive errors/reloads."""

from concurrent.futures import ThreadPoolExecutor
from threading import Event

import pytest

from omnipath_api.queries.parquet import ParquetQueries
from omnipath_api.queries.store.query_pool import QueryCapacityError, QueryPool


class Connection:
    def __init__(self):
        self.closed = False

    def close(self):
        self.closed = True


def test_capacity_reuse_failure_and_generation():
    created = []
    generation = [0]

    def factory():
        connection = Connection()
        created.append(connection)
        return connection

    pool = QueryPool(factory, lambda: generation[0], size=1, timeout=0)
    with pool.acquire() as first:
        with pytest.raises(QueryCapacityError):
            with pool.acquire():
                pytest.fail("Capacity was exceeded")
    with pytest.raises(ValueError):
        with pool.acquire() as second:
            assert second is first
            raise ValueError("Failed query")
    with pool.acquire() as third:
        assert third is first
        generation[0] += 1
    assert first.closed
    with pool.acquire() as replacement:
        assert replacement is not first
    assert len(created) == 2
    assert len(pool._connections) == 1
    pool.close()
    assert all(connection.closed for connection in created)
    assert not pool._connections


def test_pool_bounds_real_connections_and_preserves_results(tmp_path, monkeypatch):
    monkeypatch.setenv("OMNIPATH_QUERY_CONCURRENCY", "2")
    monkeypatch.setenv("OMNIPATH_QUERY_WAIT_SECONDS", "0")
    monkeypatch.setenv("OMNIPATH_DUCKDB_THREADS", "1")
    engine = ParquetQueries(tmp_path)
    started = [Event(), Event()]
    finish = Event()

    def query(index):
        with engine.query_scope():
            started[index].set()
            assert finish.wait(10)
            with engine.query_scope():  # Nested work shares its existing lease.
                return engine._db.execute("SELECT sum(i) FROM range(20) t(i)").fetchone()[0]

    try:
        with ThreadPoolExecutor(2) as executor:
            futures = [executor.submit(query, i) for i in range(2)]
            try:
                assert all(event.wait(10) for event in started)
                with pytest.raises(QueryCapacityError):
                    with engine.query_scope():
                        pytest.fail("A third database was opened")
            finally:
                finish.set()
            assert [future.result() for future in futures] == [190, 190]
        assert len(engine._query_pool._connections) == 2
        # Inventory generations retire stale databases rather than accumulating them.
        for _ in range(5):
            with engine.query_scope():
                engine.reload_resources()
                engine._reset_query_connection()
                assert engine._db.execute("SELECT 7").fetchone()[0] == 7
        assert len(engine._query_pool._connections) <= 2
    finally:
        finish.set()
        engine.close()


def test_invalid_resource_settings_fail_early(tmp_path, monkeypatch):
    monkeypatch.setenv("OMNIPATH_QUERY_CONCURRENCY", "0")
    with pytest.raises(ValueError, match="concurrency"):
        ParquetQueries(tmp_path)


def test_http_capacity_returns_retryable_status(tmp_path, monkeypatch):
    from fastapi.testclient import TestClient
    from omnipath_api.server import create_app

    monkeypatch.setenv("OMNIPATH_QUERY_CONCURRENCY", "1")
    monkeypatch.setenv("OMNIPATH_QUERY_WAIT_SECONDS", "0")
    engine = ParquetQueries(tmp_path)
    with TestClient(create_app(queries=engine)) as client:
        with engine.query_scope():
            response = client.get("/api/resources")
        assert response.status_code == 503
        assert response.headers["Retry-After"] == "1"
        assert client.get("/api/resources").status_code == 200


def test_concurrent_search_facets_and_export_match_sequential(tmp_path):
    from omnipath_core.fixtures import write_dataset

    engine = ParquetQueries(write_dataset(tmp_path))
    calls = (
        lambda: engine.search_entities_api(query="TP53", limit=10)["entities"],
        lambda: {k: v for k, v in engine.get_facets(filters={}).items() if k != "elapsed_ms"},
        lambda: engine.export_slice(filters={"genes": ["TP53"]}),
    )

    def run(call):
        with engine.query_scope(), engine.release_scope("latest"):
            return call()

    try:
        expected = [run(call) for call in calls]
        with ThreadPoolExecutor(3) as executor:
            actual = list(executor.map(run, calls))
        assert actual == expected
        assert len(engine._query_pool._connections) <= 2
    finally:
        engine.close()


def test_reload_uses_same_pool_and_never_allocates_direct_connections(tmp_path, monkeypatch):
    from omnipath_core.fixtures import write_dataset

    monkeypatch.setenv("OMNIPATH_QUERY_CONCURRENCY", "1")
    monkeypatch.setenv("OMNIPATH_QUERY_WAIT_SECONDS", "0")
    engine = ParquetQueries(write_dataset(tmp_path))
    try:
        with engine.query_scope():
            with ThreadPoolExecutor(1) as executor:
                with pytest.raises(QueryCapacityError):
                    executor.submit(engine.reload_resources).result()
        with ThreadPoolExecutor(4) as executor:
            for _ in range(4):
                executor.submit(engine.reload_resources).result()
        assert not engine._connections
        assert len(engine._query_pool._connections) <= 1
    finally:
        engine.close()


def test_busy_admin_delete_does_not_mutate_then_fail_refresh(tmp_path, monkeypatch):
    from fastapi.testclient import TestClient
    from omnipath_api.server import create_app
    from omnipath_api.settings import Settings
    from omnipath_core.fixtures import write_dataset

    monkeypatch.setenv("OMNIPATH_QUERY_CONCURRENCY", "1")
    monkeypatch.setenv("OMNIPATH_QUERY_WAIT_SECONDS", "0")
    engine = ParquetQueries(write_dataset(tmp_path))
    settings = Settings(read_only=False, admin_secret="fixture-secret")
    version = tmp_path / "resources/signor/1"
    with TestClient(create_app(queries=engine, settings=settings)) as client:
        with engine.query_scope():
            response = client.delete(
                "/api/admin/resources/signor/1", headers={"x-admin-secret": "fixture-secret"}
            )
        assert response.status_code == 503
        assert response.headers["Retry-After"] == "1"
        assert version.is_dir()
        response = client.delete(
            "/api/admin/resources/signor/1", headers={"x-admin-secret": "fixture-secret"}
        )
        assert response.status_code == 200
        assert not version.exists()
        assert not engine._connections
