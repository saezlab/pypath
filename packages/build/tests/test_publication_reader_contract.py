"""Publish bounded serving files and consume that exact manifest through PG."""

import hashlib
import json
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from omnipath_build import pipeline
from omnipath_core import PUBLISHED_TABLES, RESOURCE_FILES, SERVING_TABLES, validate_build_manifest


def test_publisher_output_is_the_shared_and_pinned_reader_contract(tmp_path, monkeypatch):
    reader = pytest.importorskip("omnipath_postgres.releases")

    def build(source, *, version, output_dir, **options):
        directory = Path(output_dir) / "resources" / source / version
        directory.mkdir(parents=True)
        result = {"resource": source, "version": version, "files": {}}
        for name, schema in {**PUBLISHED_TABLES, **SERVING_TABLES}.items():
            path = directory / f"{name}.parquet"
            pq.write_table(pa.Table.from_pylist([], schema=schema), path)
            result["files"][name] = path
        result["resolution_stats_path"] = directory / "resolution_stats.json"
        result["resolution_stats_path"].write_text("{}")
        return result

    monkeypatch.setattr(pipeline, "_build_resource", build)
    monkeypatch.setenv("__CURSOR_SANDBOX_ENV_RESTORE", "1")
    monkeypatch.setenv("HTTP_PROXY", "http://configured-proxy.invalid:8080")
    published = pipeline.build_resource(
        "signor", version="2.7", output_dir=tmp_path, max_records=20
    )
    manifest_bytes = published["manifest_path"].read_bytes()
    manifest = json.loads(manifest_bytes)
    parsed = validate_build_manifest(manifest, strict_fields=True)
    assert parsed.to_dict() == manifest
    assert manifest["max_records"] == 20
    import os

    assert os.environ["HTTP_PROXY"] == "http://configured-proxy.invalid:8080"
    # The capped build is a sample; the release must say so to be loadable.
    release_data = {
        "schema_version": 1,
        "version": "2026.10",
        "resources": {"signor": "2.7"},
        "partial_resources": True,
    }
    release_path = tmp_path / "monthly.json"
    release_path.write_text(json.dumps(release_data))
    selected = reader.load_release(tmp_path, release_path)
    (resource,) = selected.resources
    assert resource.manifest_json.encode() == manifest_bytes
    assert resource.manifest_sha256 == hashlib.sha256(manifest_bytes).hexdigest()
    assert resource.version == "2.7" and selected.version == "2026.10"
    # The reader pins the published tables; serving tables are API lookup indexes.
    assert set(resource.files) == set(RESOURCE_FILES)
    for name in RESOURCE_FILES:
        metadata = parsed.files[name]
        assert resource.files[name].rows == metadata.rows == 0
        assert resource.files[name].size_bytes == metadata.size_bytes
        assert resource.files[name].sha256 == metadata.sha256
    reader.verify_release(selected)
