"""Admin inspect/trigger/progress API tests."""

from __future__ import annotations

import json
import os
import tempfile
import time
import threading
import unittest
from unittest.mock import patch
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
from fastapi.testclient import TestClient

from omnipath_api.queries.parquet import ParquetQueries
from omnipath_api.server import create_app
from omnipath_api.jobs.worker import run_worker


class FakeOps:
    def __init__(self, delay: float = 0.0) -> None:
        self.delay = delay
        self.calls: list[tuple[str, dict]] = []

    def available(self) -> bool:
        return True

    def _maybe_wait(self, should_cancel) -> None:
        if not self.delay:
            return
        deadline = time.time() + self.delay
        while time.time() < deadline:
            if should_cancel and should_cancel():
                raise type("BuildCancelled", (Exception,), {})()
            time.sleep(0.02)

    def list_sources(self) -> list[dict]:
        return [{"source": "signor", "datasets": ["interactions"]}]

    def export_hubs(self, **kwargs):
        self.calls.append(("export_hubs", kwargs))
        on_progress = kwargs.get("on_progress")
        should_cancel = kwargs.get("should_cancel")
        if on_progress:
            on_progress(
                {
                    "pipeline": "resolver",
                    "stage": "emit:chebi",
                    "status": "running",
                    "message": "Exporting chebi",
                    "current": 0,
                }
            )
        self._maybe_wait(should_cancel)
        if on_progress:
            on_progress(
                {
                    "pipeline": "resolver",
                    "stage": "emit:chebi",
                    "status": "done",
                    "current": 12,
                    "message": "Wrote chebi (12 rows)",
                }
            )
            on_progress(
                {
                    "pipeline": "resolver",
                    "stage": "dictionaries",
                    "status": "done",
                    "message": "Wrote dictionaries",
                }
            )
            on_progress(
                {
                    "pipeline": "resolver",
                    "stage": "manifest",
                    "status": "done",
                    "message": "Wrote manifest.json",
                }
            )
        return {"chebi": 12}

    def build_all(self, **kwargs):
        self.calls.append(("build_all", kwargs))
        on_progress = kwargs.get("on_progress")
        if on_progress:
            on_progress(
                {
                    "pipeline": "resource",
                    "stage": "discover",
                    "status": "done",
                    "current": 1,
                    "message": "Found 1 dataset(s) for signor",
                }
            )
            on_progress(
                {
                    "pipeline": "resource",
                    "stage": "ingest:signor.interactions",
                    "status": "done",
                    "current": 10,
                    "message": "Read 10 records from interactions",
                }
            )
            on_progress(
                {
                    "pipeline": "resource",
                    "stage": "finalize",
                    "status": "done",
                    "message": "Completed signor/vtest",
                    "counts": {"entities": 2, "relations": 1, "payloads": 1, "rows": 10},
                }
            )
        return {"signor/vtest": {"resource": "signor", "entities_count": 2}}


class UnavailableOps:
    def available(self) -> bool:
        return False

    def list_sources(self) -> list[dict]:
        return []


def _write_lookup(path: Path) -> None:
    table = pa.table(
        {
            "source_type": ["entrez"],
            "source_id": ["7157"],
            "hub_id": ["7157"],
            "taxonomy_id": ["9606"],
            "backend": ["combined"],
        }
    )
    pq.write_table(table, path)


class TestAdminAPI(unittest.TestCase):
    def setUp(self) -> None:
        # Tests configure explicit credentials and mutation enablement.
        auth_env = patch.dict(
            os.environ,
            {
                "OMNIPATH_ADMIN_SECRET": "test-secret",
                "OMNIPATH_ADMIN_PASSWORD": "",
                "OMNIPATH_READ_ONLY": "false",
            },
        )
        auth_env.start()
        self.addCleanup(auth_env.stop)
        self.tmpdir = tempfile.TemporaryDirectory()
        self.data_root = Path(self.tmpdir.name)
        hubs = self.data_root / "reference" / "hubs"
        hubs.mkdir(parents=True)
        _write_lookup(hubs / "chebi.parquet")
        library = self.data_root / "reference" / "library"
        library.mkdir()
        (library / "manifest.json").write_text(
            json.dumps({"format": "omnipath-identity-v2", "fingerprint": "f1"})
        )
        self.ops = FakeOps()
        self.client = self._client(self.ops)

    def _client(self, ops, *, start_worker=True):
        engine = ParquetQueries(data_root=self.data_root)
        client = TestClient(create_app(queries=engine, admin_ops=ops))
        client.headers["x-admin-secret"] = "test-secret"
        self.addCleanup(engine.close)
        if start_worker:
            # Tests own and join a supervised worker; serving itself never starts one.
            self._worker_stop = threading.Event()
            self._worker = threading.Thread(
                target=run_worker,
                args=(self.data_root,),
                kwargs=dict(ops=ops, stop=self._worker_stop, poll_interval=0.01),
            )
            self._worker.start()
        return client

    def _stop_worker(self):
        self._worker_stop.set()
        self._worker.join(timeout=3)
        self.assertFalse(self._worker.is_alive())

    def tearDown(self) -> None:
        self._stop_worker()
        self.tmpdir.cleanup()

    def _wait_for_job(self, job_id: str, timeout: float = 2.0) -> dict:
        deadline = time.time() + timeout
        while time.time() < deadline:
            payload = self.client.get(f"/admin/jobs/{job_id}").json()
            if payload["status"] in {"done", "failed", "cancelled"}:
                return payload
            time.sleep(0.02)
        self.fail(f"job {job_id} did not finish")

    def test_status_inspects_resolver_files(self):
        response = self.client.get("/admin/status")
        self.assertEqual(response.status_code, 200)
        body = response.json()
        self.assertTrue(body["resolver"]["ready"])
        self.assertTrue(body["build_available"])
        self.assertEqual(body["resolver"]["library_manifest"]["fingerprint"], "f1")
        self.assertNotIn("combined", body["resolver"])

    def test_status_includes_resource_resolution_stats(self):
        version_dir = self.data_root / "resources" / "uniprot" / "v1"
        version_dir.mkdir(parents=True)
        pq.write_table(pa.table({"entity_id": [1]}), version_dir / "entity.parquet")
        pq.write_table(pa.table({"relation_id": [1]}), version_dir / "relation.parquet")
        pq.write_table(pa.table({"payload_id": [1]}), version_dir / "evidence_payloads.parquet")
        (version_dir / "resolution_stats.json").write_text(
            json.dumps(
                {
                    "scope": "unique_entity_keys",
                    "input_entities": 3,
                    "lookup_entities": 2,
                    "resolved_entities": 1,
                    "unresolved_entities": 1,
                    "not_applicable_entities": 1,
                    "by_entity_type": {},
                }
            ),
            encoding="utf-8",
        )

        built = self.client.get("/admin/status").json()["resources"]["built"]
        resource = next(item for item in built if item["resource"] == "uniprot")
        self.assertEqual(resource["resolution_stats"]["resolved_entities"], 1)
        self.assertEqual(resource["resolution_stats"]["unresolved_entities"], 1)

    def test_export_hubs_job_records_progress(self):
        started = self.client.post("/admin/jobs", json={"action": "export_hubs", "hubs": ["chebi"]})
        self.assertEqual(started.status_code, 202)
        job_id = started.json()["id"]
        self.assertTrue(any(stage["id"] == "emit:chebi" for stage in started.json()["stages"]))
        finished = self._wait_for_job(job_id)
        self.assertEqual(finished["status"], "done")
        stages = {stage["id"]: stage for stage in finished["stages"]}
        self.assertEqual(stages["emit:chebi"]["status"], "done")
        self.assertEqual(stages["emit:chebi"]["current"], 12)
        self.assertIsNotNone(stages["emit:chebi"]["started_at"])
        self.assertIsNotNone(stages["emit:chebi"]["finished_at"])
        self.assertGreaterEqual(stages["emit:chebi"]["elapsed_seconds"], 0)
        self.assertTrue(any("Wrote chebi" in line for line in finished["logs"]))
        self.assertEqual(self.ops.calls[0][0], "export_hubs")
        self.assertEqual(self.ops.calls[0][1]["hubs"], ["chebi"])
        log = self.client.get("/admin/logs/hub/chebi").json()
        self.assertTrue(log["exists"])
        self.assertIn("Wrote chebi", log["text"])
        status = self.client.get("/admin/status").json()
        chebi = next(item for item in status["resolver"]["hubs"] if item["name"] == "chebi")
        self.assertTrue(chebi["log"]["exists"])

    def test_build_resources_requires_source(self):
        response = self.client.post("/admin/jobs", json={"action": "build_resources"})
        self.assertEqual(response.status_code, 400)

    def test_build_resources_job(self):
        started = self.client.post(
            "/admin/jobs",
            json={
                "action": "build_resources",
                "sources": ["signor"],
                "version": "1.0.0",
                "max_records": 50,
            },
        )
        self.assertEqual(started.status_code, 202)
        finished = self._wait_for_job(started.json()["id"])
        self.assertEqual(finished["status"], "done")
        stages = {stage["id"]: stage for stage in finished["stages"]}
        self.assertEqual(stages["discover"]["status"], "done")
        self.assertEqual(stages["ingest:signor.interactions"]["status"], "done")
        self.assertEqual(stages["finalize"]["status"], "done")
        log = self.client.get("/admin/logs/resource/signor").json()
        self.assertTrue(log["exists"])
        self.assertIn("Completed signor/vtest", log["text"])

    def test_failed_batch_keeps_successful_resource_results(self):
        from omnipath_build.pipeline import BuildBatchError

        outcomes = {
            "signor/1": {
                "resource": "signor",
                "version": "1",
                "status": "failed",
                "error": "timeout",
            },
            "chebi/1": {"resource": "chebi", "version": "1", "entities_count": 10},
        }

        def failed_batch(**kwargs):
            raise BuildBatchError(outcomes)

        self.ops.build_all = failed_batch
        started = self.client.post(
            "/admin/jobs",
            json={
                "action": "build_resources",
                "sources": ["signor", "chebi"],
                "version": "1",
            },
        )
        self.assertEqual(started.status_code, 202)
        job = self._wait_for_job(started.json()["id"])
        self.assertEqual(job["status"], "failed")
        self.assertEqual(job["result"], outcomes)

    def test_single_flight_conflict(self):
        slow = FakeOps(delay=0.4)
        self._stop_worker()
        client = self._client(slow)
        first = client.post("/admin/jobs", json={"action": "export_hubs", "hubs": ["chebi"]})
        self.assertEqual(first.status_code, 202)
        second = client.post("/admin/jobs", json={"action": "export_hubs", "hubs": ["chebi"]})
        self.assertEqual(second.status_code, 409)
        deadline = time.time() + 2
        while time.time() < deadline:
            if client.get(f"/admin/jobs/{first.json()['id']}").json()["status"] in {
                "done",
                "failed",
                "cancelled",
            }:
                return
            time.sleep(0.02)
        self.fail("slow job did not finish")

    def test_cancel_running_job(self):
        slow = FakeOps(delay=1.0)
        self._stop_worker()
        client = self._client(slow)
        started = client.post("/admin/jobs", json={"action": "export_hubs", "hubs": ["chebi"]})
        job_id = started.json()["id"]
        cancelled = client.post(f"/admin/jobs/{job_id}/cancel")
        self.assertEqual(cancelled.status_code, 200)
        deadline = time.time() + 2
        status = None
        while time.time() < deadline:
            status = client.get(f"/admin/jobs/{job_id}").json()["status"]
            if status in {"cancelled", "done", "failed"}:
                break
            time.sleep(0.02)
        self.assertEqual(status, "cancelled")

    def test_job_events_stream(self):
        started = self.client.post("/admin/jobs", json={"action": "export_hubs", "hubs": ["chebi"]})
        job_id = started.json()["id"]
        with self.client.stream("GET", f"/admin/jobs/{job_id}/events") as stream:
            body = b"".join(stream.iter_bytes())
        self.assertIn(b"emit:chebi", body)
        self.assertIn(b"data: ", body)

    def test_admin_secret_required(self):
        previous = os.environ.get("OMNIPATH_ADMIN_SECRET")
        previous_password = os.environ.get("OMNIPATH_ADMIN_PASSWORD")
        os.environ.pop("OMNIPATH_ADMIN_PASSWORD", None)
        os.environ["OMNIPATH_ADMIN_SECRET"] = "s3cret"
        try:
            engine = ParquetQueries(data_root=self.data_root)
            client = TestClient(create_app(queries=engine, admin_ops=self.ops))
            denied = client.get("/admin/status")
            self.assertEqual(denied.status_code, 401)
            allowed = client.get("/admin/status", headers={"x-admin-secret": "s3cret"})
            self.assertEqual(allowed.status_code, 200)
        finally:
            if previous is None:
                os.environ.pop("OMNIPATH_ADMIN_SECRET", None)
            else:
                os.environ["OMNIPATH_ADMIN_SECRET"] = previous
            if previous_password is None:
                os.environ.pop("OMNIPATH_ADMIN_PASSWORD", None)
            else:
                os.environ["OMNIPATH_ADMIN_PASSWORD"] = previous_password

    def test_admin_password_env(self):
        previous = os.environ.get("OMNIPATH_ADMIN_SECRET")
        previous_password = os.environ.get("OMNIPATH_ADMIN_PASSWORD")
        os.environ.pop("OMNIPATH_ADMIN_SECRET", None)
        os.environ["OMNIPATH_ADMIN_PASSWORD"] = "hunter2"
        try:
            engine = ParquetQueries(data_root=self.data_root)
            client = TestClient(create_app(queries=engine, admin_ops=self.ops))
            self.assertEqual(client.get("/admin/status").status_code, 401)
            self.assertEqual(
                client.get("/admin/status", headers={"x-admin-secret": "hunter2"}).status_code,
                200,
            )
        finally:
            if previous is None:
                os.environ.pop("OMNIPATH_ADMIN_SECRET", None)
            else:
                os.environ["OMNIPATH_ADMIN_SECRET"] = previous
            if previous_password is None:
                os.environ.pop("OMNIPATH_ADMIN_PASSWORD", None)
            else:
                os.environ["OMNIPATH_ADMIN_PASSWORD"] = previous_password

    def test_missing_build_returns_503(self):
        client = self._client(UnavailableOps(), start_worker=False)
        response = client.post("/admin/jobs", json={"action": "export_hubs"})
        self.assertEqual(response.status_code, 503)


if __name__ == "__main__":
    unittest.main()
