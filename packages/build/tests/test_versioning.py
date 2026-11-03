"""Exercise publication and failure handling without fetching upstream datasets."""

import json
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from omnipath_build import pipeline
from omnipath_core.schema import PUBLISHED_TABLES, SERVING_TABLES


def fake_build(source, *, version, output_dir, **kwargs):
    directory = Path(output_dir) / "resources" / source / version
    directory.mkdir(parents=True)
    result = {"resource": source, "version": version, "files": {}}
    for name in (*PUBLISHED_TABLES, *SERVING_TABLES):
        pq.write_table(pa.table({"value": [1]}), directory / f"{name}.parquet")
        result["files"][name] = directory / f"{name}.parquet"
    result["resolution_stats_path"] = directory / "resolution_stats.json"
    result["resolution_stats_path"].write_text("{}")
    return result


def test_explicit_immutable_atomic_build(tmp_path, monkeypatch):
    def build(*args, **kwargs):
        assert not (tmp_path / "resources/signor/1.0.0").exists()
        return fake_build(*args, **kwargs)

    monkeypatch.setattr(pipeline, "_build_resource", build)
    result = pipeline.build_resource("signor", version="1.0.0", output_dir=tmp_path)
    assert all(path.exists() for path in result["files"].values())
    assert result["manifest_path"].exists()
    manifest = json.loads(result["manifest_path"].read_text())
    assert manifest["version"] == "1.0.0"
    assert len(manifest["files"]) == len(PUBLISHED_TABLES) + len(SERVING_TABLES)
    assert len(manifest["files"]["entity.parquet"]["sha256"]) == 64
    original = result["manifest_path"].read_bytes()
    with pytest.raises(FileExistsError):
        pipeline.build_resource("signor", version="1.0.0", output_dir=tmp_path)
    assert result["manifest_path"].read_bytes() == original
    with pytest.raises(ValueError, match="explicit numeric"):
        pipeline.build_resource("signor", output_dir=tmp_path)


def test_failed_build_is_invisible_and_retryable(tmp_path, monkeypatch):
    def fail(*args, **kwargs):
        fake_build(*args, **kwargs)
        raise RuntimeError("parser failed")

    monkeypatch.setattr(pipeline, "_build_resource", fail)
    with pytest.raises(RuntimeError):
        pipeline.build_resource("signor", version="1", output_dir=tmp_path)
    assert not (tmp_path / "resources/signor/1").exists()
    assert not list(tmp_path.glob(".build-*"))
    from omnipath_build.locking import BuildLock

    for path in (tmp_path / "resources/signor").glob("*.lock"):
        with BuildLock(path):
            pass
    monkeypatch.setattr(pipeline, "_build_resource", fake_build)
    assert pipeline.build_resource("signor", version="1", output_dir=tmp_path)["files"][
        "entity"
    ].exists()


def test_batch_pins_and_preflight(tmp_path, monkeypatch):
    monkeypatch.setattr(pipeline, "_build_resource", fake_build)
    with pytest.raises(ValueError):
        pipeline.build_all(["signor", "chebi"], versions={"signor": "1"}, output_dir=tmp_path)
    assert not (tmp_path / "resources").exists()
    result = pipeline.build_all(
        ["signor", "chebi"],
        versions={"signor": "1", "chebi": "2.1"},
        output_dir=tmp_path,
        parallel=2,
    )
    assert set(result) == {"signor/1", "chebi/2.1"}


@pytest.mark.parametrize("parallel", [1, 2])
def test_batch_continues_after_failure_and_reports_partial_results(tmp_path, monkeypatch, parallel):
    def build(source, **kwargs):
        if source == "signor":
            raise TimeoutError("upstream stopped responding")
        return fake_build(source, **kwargs)

    monkeypatch.setattr(pipeline, "_build_resource", build)
    with pytest.raises(pipeline.BuildBatchError) as error:
        pipeline.build_all(["signor", "chebi"], version="1", output_dir=tmp_path, parallel=parallel)
    assert error.value.results["signor/1"]["status"] == "failed"
    assert "TimeoutError" in error.value.results["signor/1"]["error"]
    assert (tmp_path / "resources/chebi/1/build_manifest.json").exists()
    assert not (tmp_path / "resources/signor/1").exists()
    assert not list(tmp_path.glob(".build-*"))


def test_batch_cancellation_is_not_a_skipped_resource(tmp_path, monkeypatch):
    def cancel(*args, **kwargs):
        raise pipeline.BuildCancelled()

    monkeypatch.setattr(pipeline, "build_resource", cancel)
    with pytest.raises(pipeline.BuildCancelled):
        pipeline.build_all(["signor", "chebi"], version="1", output_dir=tmp_path)
