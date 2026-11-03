"""Release pinning must hold across queries, downloads and concurrent clients."""

import hashlib
import io
import json
import os
import zipfile
from concurrent.futures import ThreadPoolExecutor

import pyarrow.parquet as pq
import pytest
from fastapi.testclient import TestClient

from omnipath_api.queries.parquet import ParquetQueries
from omnipath_api.releases import ReleaseStore
from omnipath_api.server import create_app
from omnipath_core.versioning import RESOURCE_FILES
from omnipath_core.fixtures import write_resource as write_tables


def write_resource(root, source, version, label):
    directory = root / "resources" / source / version
    directory.mkdir(parents=True)
    row = {
        key: label
        for key in (
            "relation_key",
            "subject_entity_key",
            "subject_label",
            "object_entity_key",
            "object_label",
        )
    }
    row.update(
        subject_type="protein",
        object_type="protein",
        predicate="activates",
        taxon="9606",
        category="interaction",
        interaction_class="signaling",
        is_directed=True,
        sign=1,
        sources=[source],
        evidence_count=1,
    )
    entity = dict(entity_key=label, identifier=label, label=label)
    payload = dict(relation_key=label, source=source, row_id="1", payload_json=label)
    write_tables(directory, [entity], [row], [payload])
    return directory


@pytest.fixture
def setup(tmp_path, monkeypatch):
    monkeypatch.setenv("OMNIPATH_ADMIN_SECRET", "test-release-password")
    monkeypatch.setenv("OMNIPATH_READ_ONLY", "false")
    old = write_resource(tmp_path, "signor", "1.0.0", "old")
    new = write_resource(tmp_path, "signor", "2.0.0", "new")
    write_resource(tmp_path, "extra", "1.0.0", "extra")
    store = ReleaseStore(tmp_path)
    for version, resource_version in [("2026.08", "1.0.0"), ("2026.09", "2.0.0")]:
        store.publish(
            {"schema_version": 1, "version": version, "resources": {"signor": resource_version}}
        )
    assert (tmp_path / "releases/2026.09.json").stat().st_mode & 0o777 == 0o644
    engine = ParquetQueries(tmp_path)
    with TestClient(create_app(queries=engine)) as client:
        client.headers["x-admin-secret"] = "test-release-password"
        yield store, engine, client, old, new


def test_pins_queries_and_default_independent_of_mtime(setup):
    store, engine, client, old, new = setup
    for path in old.glob("*.parquet"):
        os.utime(path, (2_000_000_000, 2_000_000_000))
    with engine.release_scope("2026.09"):
        assert engine.search_relations()["rows"][0]["subject_label"] == "new"
        with pytest.raises(ValueError):
            engine.search_relations(resources=["signor/1.0.0"])
    # Legacy default files cannot override Latest.
    (store.root / "default_release.json").write_text('{"version": "2026.08"}')
    for suffix in ("", "&release=latest", "&release=working"):
        rows = client.get("/resources?shape=svelte" + suffix).json()["resources"]
        assert {item["resource_id"]: item["version"] for item in rows} == {
            "signor": "2.0.0",
            "extra": "1.0.0",
        }
    assert client.get("/releases").json()["default"] == "latest"


def test_topic_map_follows_the_release(setup):
    _, engine, client, _, _ = setup
    for release in ("2026.08", "latest"):
        catalog = client.get(f"/resources?shape=svelte&release={release}").json()["resources"]
        expected = {}
        for resource in catalog:
            for tag in resource.get("tags") or []:
                if tag["dimension"] == "topic":
                    expected.setdefault(tag["id"], []).append(resource["resource_id"])
        topics = client.get(f"/resources/topics?release={release}").json()["topics"]
        assert {t["id"]: t["resources"] for t in topics} == {
            k: sorted(v) for k, v in expected.items()
        }
        assert all(t["label"] for t in topics)


def test_static_urls_pin_selected_resource_version(setup, monkeypatch):
    _, engine, client, _, _ = setup
    monkeypatch.setattr(engine, "public_data_url", "https://data.example.org")
    for release, version in [("2026.08", "1.0.0"), ("latest", "2.0.0")]:
        files = client.get(
            f"/resources/signor/files?release={release}&include_evidence=true"
        ).json()["files"]
        assert len(files) == len(RESOURCE_FILES)
        for file in files:
            assert (
                file["url"] == f"https://data.example.org/resources/signor/{version}/{file['name']}"
            )
        catalog = client.get(f"/resources?shape=svelte&release={release}").json()["resources"]
        signor = next(item for item in catalog if item["resource_id"] == "signor")
        assert signor["files"][0]["url"].startswith(
            f"https://data.example.org/resources/signor/{version}/"
        )
    # Existing API downloads and ZIPs remain available with static hosting enabled.
    assert client.get("/resources/signor/files/entity.parquet?release=2026.08").status_code == 200
    assert client.get("/resources/signor/download?release=2026.08").status_code == 200


def test_static_hosting_optional_and_configurable(tmp_path, monkeypatch):
    write_resource(tmp_path, "signor", "1", "test")
    monkeypatch.delenv("OMNIPATH_PUBLIC_DATA_URL", raising=False)
    engine = ParquetQueries(tmp_path)
    assert "url" not in engine.list_resource_files("signor")["files"][0]
    monkeypatch.setenv("OMNIPATH_PUBLIC_DATA_URL", "https://data.example.org/artifacts/")
    engine = ParquetQueries(tmp_path)
    assert (
        engine.list_resource_files("signor")["files"][0]["url"]
        == "https://data.example.org/artifacts/resources/signor/1/entity.parquet"
    )
    for invalid in [
        "s3://bucket",
        "https://host/?token=secret",
        "https://user:pass@host",
        "https://host/#fragment",
    ]:
        with pytest.raises(ValueError, match="OMNIPATH_PUBLIC_DATA_URL"):
            ParquetQueries(tmp_path, public_data_url=invalid)


def test_latest_numeric_order_and_new_resources(setup):
    store, engine, client, _, _ = setup
    # A later-written lower version must not win; numeric ordering beats lexical ordering.
    write_resource(store.root, "extra", "1.10.0", "ten")
    write_resource(store.root, "extra", "1.9.0", "nine")
    legacy = write_resource(store.root, "extra", "abcdef123456", "legacy")
    for path in legacy.glob("*.parquet"):
        os.utime(path, (2_000_000_000, 2_000_000_000))
    write_resource(store.root, "new_source", "1", "added")
    pending = write_resource(store.root, "signor", "3.0.0", "pending")
    table = pending / "relation_evidence.parquet"
    table.rename(pending / "relation_evidence.partial")
    rows = client.get("/resources?shape=svelte").json()["resources"]
    assert {item["resource_id"]: item["version"] for item in rows} == {
        "signor": "2.0.0",
        "extra": "1.10.0",
        "new_source": "1",
    }
    (pending / "relation_evidence.partial").rename(table)
    with engine.release_scope("latest"):
        assert engine.resolve_resource_info("signor")["version"] == "3.0.0"
    with engine.release_scope("2026.09"):
        assert engine.resolve_resource_info("signor")["version"] == "2.0.0"


def test_downloads_and_parallel_request_isolation(setup):
    _, _, client, _, _ = setup

    def read(version):
        response = client.get(
            "/api/resources/signor/files/entity.parquet", headers={"X-OmniPath-Release": version}
        )
        assert response.status_code == 200
        return pq.read_table(io.BytesIO(response.content))["label"][0].as_py()

    with ThreadPoolExecutor(4) as pool:
        assert list(pool.map(read, ["2026.08", "2026.09"] * 5)) == ["old", "new"] * 5
    assert client.get("/resources/extra/download?release=2026.09").status_code == 404
    assert client.get("/resources?release=missing").status_code == 404
    assert client.get("/resources?release=../../escape").status_code == 400


def _zip_names(client, url, **kwargs):
    response = client.get(url, **kwargs)
    assert response.status_code == 200
    with zipfile.ZipFile(io.BytesIO(response.content)) as archive:
        return set(archive.namelist())


def test_downloads_omit_evidence_by_default(setup):
    _, _, client, _, _ = setup
    catalog = client.get("/resources?shape=svelte").json()["resources"]
    signor = next(item for item in catalog if item["resource_id"] == "signor")
    published = set(RESOURCE_FILES) - {"evidence_payloads.parquet"}
    assert {item["name"] for item in signor["files"]} == published

    listing = client.get("/resources/signor/files").json()["files"]
    assert {item["name"] for item in listing} == published
    with_evidence = client.get("/resources/signor/files?include_evidence=true").json()["files"]
    assert {item["name"] for item in with_evidence} == set(RESOURCE_FILES)

    assert _zip_names(client, "/resources/signor/download") == published
    assert _zip_names(client, "/resources/signor/download?include_evidence=true") == set(
        RESOURCE_FILES
    )

    bundled = client.post("/resources/download", json={"resource_ids": ["signor"]})
    assert bundled.status_code == 200
    with zipfile.ZipFile(io.BytesIO(bundled.content)) as archive:
        assert set(archive.namelist()) == published

    evidence = client.get("/resources/signor/files/evidence_payloads.parquet")
    assert evidence.status_code == 200
    assert evidence.content.startswith(b"PAR1")


def test_invalid_releases_and_pinned_deletion(setup):
    store, _, client, old, _ = setup
    with pytest.raises(FileExistsError):
        store.publish(store.get("2026.08"))
    with pytest.raises(ValueError, match="Missing release artifact"):
        store.publish({"schema_version": 1, "version": "2027", "resources": {"missing": "1"}})
    assert not (store.root / "releases" / "2027.json").exists()
    response = client.delete("/admin/resources/signor/versions/1.0.0")
    assert response.status_code == 400
    assert "pinned" in response.text
    assert old.exists()
    (old / "relation.parquet").unlink()
    response = client.get("/resources?release=2026.08")
    assert response.status_code == 400  # Never fall back to the newer resource.


def test_admin_publish_and_explicit_build_version(setup):
    _, _, client, _, _ = setup
    payload = {"schema_version": 1, "version": "2027", "resources": {"signor": "2.0.0"}}
    assert client.post("/admin/releases", json=payload).status_code == 201
    assert client.get("/releases").json()["default"] == "latest"
    assert client.post("/admin/releases", json=payload).status_code == 409
    assert (
        client.post(
            "/admin/jobs", json={"action": "build_resource", "source": "signor"}
        ).status_code
        == 400
    )
    assert (
        client.post(
            "/admin/jobs", json={"action": "build_resource", "source": "signor", "version": "1.0.0"}
        ).status_code
        == 409
    )


def _with_manifest(root, source, version, **fields):
    write_resource(root, source, version, source)
    path = root / "resources" / source / version / "build_manifest.json"
    path.write_text(json.dumps({"resource": source, "version": version, **fields}))
    return path


def test_publish_pins_build_manifest_bytes(tmp_path):
    path = _with_manifest(tmp_path, "signor", "1", max_records=None)
    store = ReleaseStore(tmp_path)
    with pytest.warns(UserWarning, match="taxonomy"):
        published = store.publish(
            {"schema_version": 1, "version": "1", "resources": {"signor": "1"}}
        )
    assert published["resource_manifests"] == {
        "signor": hashlib.sha256(path.read_bytes()).hexdigest()
    }
    with pytest.raises(ValueError, match="do not match"):
        store.publish(
            {
                "schema_version": 1,
                "version": "2",
                "resources": {"signor": "1"},
                "resource_manifests": {"signor": "0" * 64},
            }
        )


@pytest.mark.parametrize("fields", [{"max_records": 20}, {"datasets": ["interactions"]}])
def test_publish_refuses_unmarked_sample_builds(tmp_path, fields):
    _with_manifest(tmp_path, "signor", "1", **fields)
    store = ReleaseStore(tmp_path)
    release = {"schema_version": 1, "version": "1", "resources": {"signor": "1"}}
    with pytest.raises(ValueError, match="sample builds"):
        store.publish(release)
    with pytest.warns(UserWarning, match="taxonomy"):
        assert store.publish({**release, "partial_resources": True})["partial_resources"] is True
