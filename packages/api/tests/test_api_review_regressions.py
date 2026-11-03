"""Behavioral regression cases from the September 2026 code review."""

from __future__ import annotations

import io
import json
import sqlite3
import threading
from concurrent.futures import ThreadPoolExecutor

import pyarrow.parquet as pq
import pytest
from fastapi.testclient import TestClient

from omnipath_api.admin import AdminService
from omnipath_api.queries.parquet import ParquetQueries
from omnipath_api.jobs.runner import Job
from omnipath_api.jobs.store import JobStore
from omnipath_api.jobs.worker import run_worker, worker_lock
from omnipath_api.server import create_app
from omnipath_api.settings import Settings
from omnipath_core.fixtures import fixture_key as key, write_dataset
from omnipath_core.fixtures import write_resource


@pytest.fixture
def engine(tmp_path):
    result = ParquetQueries(write_dataset(tmp_path))
    yield result
    result.close()


@pytest.mark.parametrize(
    "field",
    ["entity_ids", "entityIds", "entities", "genes", "entity_pks", "entityPks", "scope_entity_ids"],
)
def test_zero_match_filters_never_broaden_relations(engine, field):
    filters = {field: ["uniprot|MISSING"]}
    assert engine.search_relations(filters)["total"] == 0
    assert engine.get_facets(filters)["categories"] == {}
    assert engine.get_scoped_relation_facets(filters) == []
    exported, _ = engine.export_slice(filters)
    assert pq.read_table(io.BytesIO(exported)).num_rows == 0


@pytest.mark.parametrize("token", ["P04637", "uniprot|P04637", "TP53"])
def test_short_accession_label_and_public_id_resolve_consistently(engine, token):
    filters = {"entityIds": [token]}
    assert engine.search_relations(filters)["total"] == 1
    assert engine.get_facets(filters)["categories"] == {"interaction": 1}
    assert any(
        item["facetValue"] == "affects" for item in engine.get_scoped_relation_facets(filters)
    )
    exported, _ = engine.export_slice(filters)
    assert pq.read_table(io.BytesIO(exported))["relation_key"].to_pylist() == [key("relation")]


@pytest.mark.parametrize(
    "query,phase",
    [("", "all"), ("Match", "prefix"), ("Inside", "contains"), ("ALIAS", "prefix")],
)
def test_entity_cursor_traverses_complete_sort_order(tmp_path, query, phase):
    directory = tmp_path / "resources" / "test" / "1"
    directory.mkdir(parents=True)
    rows = [
        dict(
            entity_key=f"{100 - i:064x}",
            namespace="uniprot",
            identifier=f"P{i}",
            label=label,
            entity_type="protein",
            taxon="9606",
            has_hierarchy=False,
            parent_count=0,
            child_count=0,
            identifiers=[dict(ns="name", id=f"ALIAS{i}", source="test")],
            annotations=[],
        )
        for i, label in enumerate(
            ["Match", "Match A", "Match B", "X Inside A", "X Inside B", "X Inside C"]
        )
    ]
    write_resource(directory, rows, [])
    engine = ParquetQueries(tmp_path)
    try:
        expected = engine.search_entities_api(query, limit=100)["entities"]
        cursor = None
        found = []
        for _ in range(20):
            page = engine.search_entities_api(query, limit=1, cursor=cursor)
            found.extend(page["entities"])
            cursor = page["nextCursor"]
            if cursor is None:
                break
            assert cursor["phase"] == phase
        assert [e["entityPk"] for e in found] == [e["entityPk"] for e in expected]
        assert len({e["entityPk"] for e in found}) == len(found)
        first_cursor = engine.search_entities_api(query, limit=1)["nextCursor"]
        with pytest.raises(ValueError, match="another query"):
            engine.search_entities_api(query + " changed", cursor=first_cursor)
    finally:
        engine.close()


def test_api_rejects_malformed_filters_and_unicode_alias_paths(engine):
    with TestClient(create_app(queries=engine)) as client:
        for raw in ["{", "[]", "42", "null", '{"entity_ids":"P04637"}', '{"gene_mode":"invalid"}']:
            assert client.get("/api/entities/search", params={"filters": raw}).status_code == 422
        assert (
            client.post("/relations/search", json={"filters": {"entity_ids": "bad"}}).status_code
            == 422
        )
        assert client.get("/api/entities/λ").status_code == 404
        assert client.get("/app-api/entities/café").status_code == 404
        first = client.get("/api/entities/search", params={"limit": 1}).json()
        assert first["nextCursor"] is not None
        assert (
            client.get(
                "/api/entities/search",
                params={"limit": 1, "cursor": json.dumps(first["nextCursor"])},
            ).status_code
            == 200
        )


def test_admin_is_disabled_without_secret_and_readonly_blocks_every_mutation(engine):
    client = TestClient(
        create_app(queries=engine, settings=Settings(admin_secret="", read_only=False))
    )
    assert client.get("/admin/status").status_code == 503
    assert client.post("/admin/jobs", json={"action": "export_hubs"}).status_code == 503
    client = TestClient(
        create_app(queries=engine, settings=Settings(admin_secret="secret", read_only=True))
    )
    assert client.get("/admin/status").status_code == 401
    client.headers["x-admin-secret"] = "secret"
    for path, body in [
        ("/admin/jobs", {"action": "export_hubs"}),
        ("/admin/releases", {}),
        ("/admin/jobs/abc/cancel", None),
    ]:
        assert client.post(path, json=body).status_code == 403
    assert client.delete("/admin/resources/test/1").status_code == 403
    assert (
        client.get("/admin/status?secret=secret", headers={"x-admin-secret": ""}).status_code == 401
    )


def test_inventory_generation_stays_pinned_through_reload(engine, tmp_path):
    old = engine.inventory_snapshot()
    with engine.release_scope("latest", snapshot=old):
        before = engine._table_paths("entity")
        target = tmp_path / "resources" / "uniprot" / "2"
        import shutil

        shutil.copytree(target.parent / "1", target)
        with ThreadPoolExecutor(1) as pool:
            pool.submit(engine.reload_resources).result()
        assert engine._table_paths("entity") == before
        assert "uniprot/2" not in engine.resources
    assert "uniprot/2" in engine.resources
    with pytest.raises(TypeError):
        old.resources["uniprot/1"]["version"] = "oops"


def test_request_reuses_snapshot_for_multiple_engine_calls(engine, tmp_path):
    from fastapi import Request
    from omnipath_api.routers.common import _call

    app = create_app(queries=engine)
    request = Request(
        {
            "type": "http",
            "method": "GET",
            "path": "/",
            "headers": [],
            "query_string": b"",
            "app": app,
        }
    )
    before = _call(request, "list_resources")
    import shutil

    shutil.copytree(tmp_path / "resources/uniprot/1", tmp_path / "resources/uniprot/2")
    engine.reload_resources()
    assert _call(request, "list_resources") == before
    assert len(engine.list_resources()) == len(before) + 1


def test_ontology_scope_ambiguity_and_all_child_pages(tmp_path):
    target = tmp_path / "resources/test/1"
    rows = [
        dict(
            entity_key=f"{namespace}:{i}",
            namespace=namespace,
            identifier=f"T{i}",
            label=f"Term {i}",
        )
        for namespace in ["one", "two"]
        for i in range(106)
    ]
    relations = [
        dict(
            relation_key=f"{namespace}:{i}:0",
            subject_entity_key=f"{namespace}:{i}",
            subject_label=f"Term {i}",
            object_entity_key=f"{namespace}:0",
            object_label="Term 0",
            predicate="subclass_of",
            statement_kind="ontology",
        )
        for namespace in ["one", "two"]
        for i in range(1, 106)
    ]
    relations.append(
        dict(
            relation_key="two:1:one:0",
            subject_entity_key="two:1",
            subject_label="Other",
            object_entity_key="one:0",
            object_label="Term 0",
            predicate="subclass_of",
            statement_kind="ontology",
        )
    )
    write_resource(target, rows, relations)
    engine = ParquetQueries(tmp_path)
    try:
        with pytest.raises(ValueError, match="Ambiguous"):
            engine.get_ontology_children("T0")
        seen = []
        offset = 0
        while True:
            result = engine.get_ontology_children("T0", ontology_id="one", offset=offset, limit=20)
            assert result["total"] == 105
            seen.extend(result["children"])
            if result["nextCursor"] is None:
                break
            offset = result["nextCursor"]
        assert len({r["entityPk"] for r in seen}) == 105
        assert all(r["ontologyId"] == "one" for r in seen)
        assert engine.get_ontology_tree(["one|T0"])["root"]["hasMoreChildren"]
        assert engine.get_ontology_children("one|T0")["total"] == 105
        with pytest.raises(ValueError, match="contradicts"):
            engine.get_ontology_tree(["one|T0"], ontology_id="two")
    finally:
        engine.close()


class WorkerOps:
    def available(self):
        return True

    def export_hubs(self, **kwargs):
        return {"completed": True}


def test_durable_jobs_survive_api_restart_and_worker_is_separate(tmp_path):
    first = AdminService(tmp_path, ops=WorkerOps())
    initial = first.start("export_hubs")
    assert initial["status"] == "queued"
    assert not any(t.name.startswith("admin-") for t in threading.enumerate())
    restarted = AdminService(tmp_path, ops=WorkerOps())
    assert restarted.get_job(initial["id"])["status"] == "queued"
    with pytest.raises(RuntimeError, match="already"):
        restarted.start("export_hubs")
    run_worker(tmp_path, ops=WorkerOps(), once=True)
    assert restarted.get_job(initial["id"])["status"] == "done"


def test_cancel_after_completed_work_still_records_result(tmp_path):
    service = AdminService(tmp_path, ops=WorkerOps())
    job = Job("export_hubs", {}, [])
    service.run_job(job, should_stop=lambda: True)
    assert job.status == "done"
    assert job.result == {"completed": True}


def test_worker_lock_recovery_and_cancelled_queue(tmp_path):
    store = JobStore(tmp_path)
    job = Job("export_hubs", {}, [])
    store.enqueue(job)
    claimed = store.claim()
    assert claimed.id == job.id
    with worker_lock(tmp_path):
        with pytest.raises(RuntimeError, match="Another"):
            run_worker(tmp_path, ops=WorkerOps(), once=True)
    run_worker(tmp_path, ops=WorkerOps(), once=True)
    assert store.get(job.id).status == "failed"
    assert "stopped before completion" in store.get(job.id).error
    next_job = Job("export_hubs", {}, [])
    store.enqueue(next_job)
    assert store.cancel(next_job.id).status == "cancelled"
    assert store.claim() is None


def test_queue_persistence_failure_is_not_success(tmp_path, monkeypatch):
    service = AdminService(tmp_path, ops=WorkerOps())

    def fail(_job):
        raise OSError("disk full")

    monkeypatch.setattr(service.job_store, "enqueue", fail)
    with pytest.raises(OSError, match="disk full"):
        service.start("export_hubs")


def test_queue_readers_tolerate_file_created_before_its_table(tmp_path):
    store = JobStore(tmp_path)
    store.path.parent.mkdir(parents=True)
    # The writer's connect() creates the file before CREATE TABLE runs.
    sqlite3.connect(store.path).close()
    assert store.list() == []
    assert store.get("missing") is None


def test_readonly_startup_and_queries_do_not_create_files(tmp_path):
    root = write_dataset(tmp_path)
    before = set(root.rglob("*"))
    with TestClient(create_app(data_root=root, settings=Settings(read_only=True))) as client:
        assert client.get("/entities/search?q=ADA").status_code == 200
        assert client.post("/relations/search", json={}).status_code == 200
        assert client.get("/entities/examples").status_code == 200
    assert set(root.rglob("*")) == before


def test_nested_facet_filters_preserve_outer_entity_selection(engine):
    assert (
        engine.get_scoped_relation_facets(
            {"entityIds": ["uniprot|MISSING"], "filters": {"predicates": ["affects"]}}
        )
        == []
    )
