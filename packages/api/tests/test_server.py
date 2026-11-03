"""Integration tests for the FastAPI parquet serving handlers."""

import io
import unittest
import tempfile
import zipfile
from omnipath_core.fixtures import write_dataset
from omnipath_api.settings import Settings

from fastapi.testclient import TestClient

from omnipath_api.queries.parquet import ParquetQueries
from omnipath_api.server import create_app


class TestServerAPI(unittest.TestCase):
    client: TestClient

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.tmp.cleanup)
        engine = ParquetQueries(data_root=write_dataset(cls.tmp.name))
        cls.addClassCleanup(engine.close)
        cls.client = TestClient(create_app(queries=engine, settings=Settings(admin_secret="test")))
        cls.client.headers["x-admin-secret"] = "test"

    def test_admin_status(self):
        # Handler-shape test: independent of developer .env authentication settings.
        response = self.client.get("/admin/status")
        self.assertEqual(response.status_code, 200)
        body = response.json()
        self.assertIn("resolver", body)
        self.assertIn("resources", body)
        self.assertIn("job", body)

    def test_health(self):
        response = self.client.get("/health")
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json().get("status"), "ok")

    def test_health_api_prefix(self):
        response = self.client.get("/api/health")
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json().get("status"), "ok")

    def test_resources_catalog(self):
        response = self.client.get("/resources/catalog")
        self.assertEqual(response.status_code, 200)
        self.assertIn("resources", response.json())

    def test_entities_search(self):
        response = self.client.get("/entities/search?q=TP53")
        self.assertEqual(response.status_code, 200)
        self.assertIn("entities", response.json())

    def test_relations_search(self):
        response = self.client.post("/relations/search", json={"filters": {}, "limit": 5})
        self.assertEqual(response.status_code, 200)
        data = response.json()
        self.assertIn("relations", data)
        self.assertIn("total", data)

    def test_openapi_document(self):
        response = self.client.get("/openapi.json")
        self.assertEqual(response.status_code, 200)
        spec = response.json()
        self.assertEqual(spec["info"]["title"], "OmniPath Parquet API")
        paths = spec["paths"]
        self.assertIn("/health", paths)
        self.assertIn("/entities/search", paths)
        self.assertIn("/relations/search", paths)
        self.assertIn("/resources", paths)
        self.assertIn("/admin/status", paths)
        self.assertIn("/admin/jobs", paths)
        servers = spec.get("servers") or []
        self.assertTrue(any("/api" in str(server.get("url") or "") for server in servers))

    def test_docs_page(self):
        response = self.client.get("/docs")
        self.assertEqual(response.status_code, 200)
        self.assertIn("text/html", response.headers["content-type"])
        self.assertIn("swagger", response.text.lower())

    def test_docs_page_api_prefix(self):
        response = self.client.get("/api/docs")
        self.assertEqual(response.status_code, 200)
        self.assertIn("text/html", response.headers["content-type"])

    def test_prefixed_app_api_alias(self):
        response = self.client.get("/app-api/health")
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json().get("status"), "ok")

    def test_resource_files_and_downloads(self):
        catalog = self.client.get("/resources?shape=svelte").json()["resources"]
        self.assertGreater(len(catalog), 0)
        resource_id = catalog[0]["resource_id"]
        self.assertTrue(catalog[0].get("files"))

        catalog_names = {item["name"] for item in catalog[0]["files"]}
        self.assertIn("entity.parquet", catalog_names)
        self.assertNotIn("evidence_payloads.parquet", catalog_names)

        listing = self.client.get(f"/resources/{resource_id}/files")
        self.assertEqual(listing.status_code, 200)
        files = listing.json()["files"]
        names = {item["name"] for item in files}
        self.assertIn("entity.parquet", names)
        self.assertNotIn("evidence_payloads.parquet", names)
        self.assertTrue(any(item.get("columns") for item in files))

        parquet = self.client.get(f"/resources/{resource_id}/files/entity.parquet")
        self.assertEqual(parquet.status_code, 200)
        self.assertTrue(parquet.content.startswith(b"PAR1"))

        missing = self.client.get(f"/resources/{resource_id}/files/../engine.py")
        self.assertEqual(missing.status_code, 404)

        zipped = self.client.get(f"/resources/{resource_id}/download")
        self.assertEqual(zipped.status_code, 200)
        self.assertIn("zip", zipped.headers.get("content-type", ""))
        self.assertGreater(len(zipped.content), 4)
        self.assertEqual(zipped.content[:2], b"PK")
        with zipfile.ZipFile(io.BytesIO(zipped.content)) as archive:
            self.assertNotIn("evidence_payloads.parquet", archive.namelist())


if __name__ == "__main__":
    unittest.main()
