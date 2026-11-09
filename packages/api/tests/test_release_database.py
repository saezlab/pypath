"""Serving releases from their databases: build, pick up, prune, and what a pruned
release answers."""

import json

import duckdb
import pytest
from fastapi.testclient import TestClient

from omnipath_api.cli import main as cli
from omnipath_api.queries.release import ReleaseQueries
from omnipath_api.release_database import (
    build_served_database,
    prune_served_databases,
    publish_newest,
    served_databases,
)
from omnipath_api.server import create_app
from omnipath_core.fixtures import write_dataset


@pytest.fixture
def root(tmp_path):
    root = write_dataset(tmp_path / "data")
    for version in ("1", "2", "3"):
        publish_newest(root, version)
    return root


def test_build_records_the_database_and_prunes_old_ones(root):
    for version in ("1", "2"):
        build_served_database(root, version, keep=None, memory_limit="1GB", threads=1)
    record = build_served_database(root, "3", keep=2, memory_limit="1GB", threads=1)
    assert record["pruned"] == ["1"]
    assert served_databases(root) == ["3", "2"]
    written = json.loads((root / "serving" / "3.json").read_text())
    assert written["release"] == "3" and written["resources"] == 3
    assert written["size_bytes"] == (root / "serving" / "3.duckdb").stat().st_size
    assert not (root / "serving" / "1.json").exists()
    # The release manifests stay: older releases remain downloadable.
    assert (root / "releases" / "1.json").is_file()
    with duckdb.connect(str(root / "serving" / "3.duckdb"), read_only=True) as db:
        assert db.execute("SELECT release FROM database_info").fetchone() == ("3",)
    assert prune_served_databases(root, 1) == ["2"]
    with pytest.raises(ValueError):
        prune_served_databases(root, 0)


def test_unknown_release_is_refused(root):
    with pytest.raises(FileNotFoundError):
        build_served_database(root, "9", memory_limit="1GB", threads=1)
    assert served_databases(root) == []


def test_latest_is_the_newest_served_release_and_pruned_ones_answer_410(root):
    engine = ReleaseQueries(root)
    with TestClient(create_app(queries=engine)) as client:
        # Nothing served yet.
        assert client.post("/relations/search", json={}).status_code == 503
        build_served_database(root, "1", keep=None, memory_limit="1GB", threads=1)
        assert client.post("/relations/search", json={}).json()["total"] == 1
        build_served_database(root, "2", keep=1, memory_limit="1GB", threads=1)
        with engine.release_scope("latest") as version:
            assert version == "2"
        response = client.post("/relations/search?release=1", json={})
        assert response.status_code == 410
        assert "no longer served" in response.json()["detail"]
        assert client.post("/relations/search?release=2", json={}).status_code == 200
        # A release that was never published.
        assert client.post("/relations/search?release=7", json={}).status_code == 404
        # A published release without a database is not served either.
        assert client.post("/relations/search?release=3", json={}).status_code == 410
    engine.close()


def test_cli_builds_and_prunes(root, capsys):
    assert cli(["build-release-db", "1", "--data-root", str(root), "--threads", "1", "--memory-limit", "1GB"]) == 0
    assert json.loads(capsys.readouterr().out)["release"] == "1"
    assert cli(["build-release-db", "2", "--data-root", str(root), "--threads", "1", "--memory-limit", "1GB", "--keep", "0"]) == 0
    capsys.readouterr()
    assert served_databases(root) == ["2", "1"]
    assert cli(["prune-release-dbs", "--data-root", str(root), "--keep", "1"]) == 0
    assert json.loads(capsys.readouterr().out) == {"removed": ["1"]}
