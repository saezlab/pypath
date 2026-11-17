"""Publishing a release builds its subsets before the database it is served from."""

import json

import pytest

from omnipath_api import publish_release, release_database
from omnipath_core.fixtures import write_dataset


@pytest.fixture
def manifest(tmp_path):
    root = write_dataset(tmp_path / "data")
    resources = {path.parent.name: path.name for path in (root / "resources").glob("*/*")}
    path = tmp_path / "release.json"
    path.write_text(json.dumps({"schema_version": 1, "version": "5", "resources": resources}))
    return root, path


@pytest.fixture
def calls(monkeypatch):
    seen = []
    monkeypatch.setattr(
        "omnipath_subsets.parquet_products.build_release_subsets",
        lambda root, version, *args, **kwargs: seen.append(("subsets", version)) or {},
    )
    monkeypatch.setattr(
        release_database,
        "build_served_database",
        lambda root, version, **kwargs: seen.append(("database", version)) or {},
    )
    return seen


def test_the_subsets_are_built_before_the_serving_database(manifest, calls):
    root, path = manifest
    publish_release.main([str(path), "--data-root", str(root), "--build-database"])
    assert calls == [("subsets", "5"), ("database", "5")]


def test_subsets_can_be_skipped_explicitly(manifest, calls):
    root, path = manifest
    publish_release.main(
        [str(path), "--data-root", str(root), "--build-database", "--skip-subsets"]
    )
    assert calls == [("database", "5")]


def test_a_release_published_without_a_database_builds_no_subsets(manifest, calls):
    root, path = manifest
    publish_release.main([str(path), "--data-root", str(root)])
    assert calls == []


def test_a_failing_subset_build_stops_before_the_database(manifest, calls, monkeypatch):
    root, path = manifest

    def fail(*args, **kwargs):
        raise RuntimeError("subsets failed")

    monkeypatch.setattr("omnipath_subsets.parquet_products.build_release_subsets", fail)
    with pytest.raises(RuntimeError, match="subsets failed"):
        publish_release.main([str(path), "--data-root", str(root), "--build-database"])
    assert calls == []
