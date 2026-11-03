"""The serving status: query slots, recent response times, coarse load."""

import threading

from fastapi.testclient import TestClient

from omnipath_api.queries.parquet import ParquetQueries
from omnipath_api.server import create_app
from omnipath_api.status import RequestTimes
from omnipath_core.fixtures import write_resource


def test_status_reports_slots_and_recent_requests(tmp_path):
    write_resource(tmp_path / "resources/a/1", [dict(entity_key="p1", label="TP53")])
    client = TestClient(create_app(queries=ParquetQueries(tmp_path)))
    first = client.get("/status").json()
    assert first["state"] == "ok"
    assert first["queries"]["running"] == 0 and first["queries"]["waiting"] == 0
    assert first["queries"]["slots"] >= 1
    assert first["requestsLastMinute"] == 0  # status and health requests are not timed
    client.get("/entities/search?q=TP53")
    second = client.get("/status").json()
    assert second["requestsLastMinute"] == 1 and second["medianResponseMs"] >= 0


def test_waiting_queries_mark_the_service_queued(tmp_path):
    write_resource(tmp_path / "resources/a/1", [dict(entity_key="p1", label="TP53")])
    engine = ParquetQueries(tmp_path)
    pool = engine._query_pool
    from omnipath_api.status import status

    held, release = [], threading.Event()

    def hold():
        with pool.acquire():
            held.append(True)
            release.wait(5)

    workers = [threading.Thread(target=hold) for _ in range(pool.size + 1)]
    for worker in workers:
        worker.start()
    try:
        for _ in range(200):
            if pool.waiting:
                break
            threading.Event().wait(0.01)
        current = status(pool, RequestTimes())
        assert current["state"] == "queued"
        assert current["queries"] == {"running": pool.size, "slots": pool.size, "waiting": 1}
    finally:
        release.set()
        for worker in workers:
            worker.join()
    assert (pool.running, pool.waiting) == (0, 0)
