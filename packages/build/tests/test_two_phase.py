import importlib
import multiprocessing

import duckdb
import pytest

from omnipath_build.discovery import DiscoveredDataset
from omnipath_build.pipeline import build_resource
from omnipath_build.progress import BuildCancelled


@pytest.fixture
def dataset(tmp_path, monkeypatch):
    module = tmp_path / "batch_fixture_source.py"
    module.write_text("""
import time
def mapper(row):
    if row.get('fail'):
        raise ValueError('deliberate mapper failure')
    if row.get('slow'):
        time.sleep(20)
    if row.get('skip'):
        return None
    return {'subject': {'type': 'protein', 'identifiers': [
        {'type':'uniprot','value': 'P' + str(row['i'] % 3)},
        {'type':'genesymbol','value': 'alias' + str(row['i'] % 7)}]},
        'predicate': 'affects',
        'object': {'type':'protein','identifiers':[{'type':'uniprot','value':'P9'}]},
        'annotations':[{'term':'description','value':str(row['i'] % 5)}]}
class Input:
    mapper = staticmethod(mapper)
    def raw(self, **kwargs):
        row = {}
        for i in range(123):
            row['i'] = i
            row['skip'] = i % 17 == 0
            yield row
items = Input()
""")
    monkeypatch.syspath_prepend(str(tmp_path))
    import sys

    sys.modules.pop("batch_fixture_source", None)
    mod = importlib.import_module("batch_fixture_source")
    ds = DiscoveredDataset(
        "fixture", "items", mod.__name__, lambda: (), raw_dataset=mod.items, mapper=mod.mapper
    )
    monkeypatch.setattr(
        "omnipath_build.pipeline.discover_datasets", lambda **kw: ("fixture", [ds], {})
    )
    return ds


def test_worker_counts_preserve_outputs_statistics_and_row_limits(tmp_path, dataset):
    empty = tmp_path / "library"
    results = []
    for workers in (1, 2):
        results.append(
            build_resource(
                "fixture",
                version="1",
                output_dir=tmp_path / str(workers),
                batch_workers=workers,
                max_records=101,
                max_batch_records=19,
                batch_size=11,
                max_batch_relations=9,
                max_batch_bytes=10000,
                library_dir=empty,
                duckdb_memory_limit="256MB",
                progress=False,
            )
        )
    for key in (
        "input_entities",
        "resolved_entities",
        "unresolved_entities",
        "by_rule",
        "by_entity_type",
    ):
        assert results[0]["resolution_stats"][key] == results[1]["resolution_stats"][key]
    with duckdb.connect() as db:
        for name in results[0]["files"]:
            a, b = (str(r["files"][name]) for r in results)
            assert (
                db.execute(
                    """SELECT count(*) FROM (
                (SELECT * FROM read_parquet(?) EXCEPT ALL SELECT * FROM read_parquet(?)) UNION ALL
                (SELECT * FROM read_parquet(?) EXCEPT ALL SELECT * FROM read_parquet(?)))""",
                    [a, b, b, a],
                ).fetchone()[0]
                == 0
            )
        assert (
            db.execute(
                "SELECT count(*) FROM read_parquet(?)",
                [str(results[1]["files"]["evidence_payloads"])],
            ).fetchone()[0]
            == 95
        )
    assert not list(tmp_path.rglob("phase-*"))
    assert not multiprocessing.active_children()


def test_two_phase_worker_failure_does_not_publish(tmp_path, dataset, monkeypatch):
    monkeypatch.setattr(dataset.raw_dataset, "raw", lambda **kwargs: iter([{"i": 1, "fail": True}]))
    with pytest.raises(
        RuntimeError, match=r"fixture:items input row 0 .*deliberate mapper failure"
    ):
        build_resource(
            "fixture",
            version="1",
            output_dir=tmp_path / "failed",
            batch_workers=2,
            library_dir=tmp_path / "empty",
            progress=False,
        )
    assert not (tmp_path / "failed/resources/fixture/1").exists()
    assert not list((tmp_path / "failed").glob(".build-*"))
    assert not multiprocessing.active_children()


def test_dataset_without_rows_fails_unless_it_may_be_empty(tmp_path, dataset, monkeypatch):
    monkeypatch.setattr(dataset.raw_dataset, "raw", lambda **kwargs: iter(()))
    build = lambda name: build_resource(
        "fixture",
        version="1",
        output_dir=tmp_path / name,
        batch_workers=1,
        library_dir=tmp_path / "empty",
        progress=False,
    )
    with pytest.raises(RuntimeError, match=r"fixture:items yielded no rows"):
        build("failed")
    assert not (tmp_path / "failed/resources/fixture/1").exists()
    monkeypatch.setattr(dataset.raw_dataset, "allow_empty", True, raising=False)
    build("allowed")
    assert (tmp_path / "allowed/resources/fixture/1").exists()
    assert not multiprocessing.active_children()


def test_two_phase_cancellation_reaps_workers_and_preserves_atomicity(tmp_path, dataset):
    cancelled = False

    def progress(event):
        nonlocal cancelled
        if event.get("stage") == "prepare":
            cancelled = True

    with pytest.raises(BuildCancelled):
        build_resource(
            "fixture",
            version="1",
            output_dir=tmp_path / "cancelled",
            batch_workers=2,
            library_dir=tmp_path / "empty",
            progress=False,
            on_progress=progress,
            should_cancel=lambda: cancelled,
        )
    assert not (tmp_path / "cancelled/resources/fixture/1").exists()
    assert not multiprocessing.active_children()


def test_default_build_publishes_two_phase_manifest(tmp_path, dataset):
    import json

    result = build_resource(
        "fixture",
        version="1",
        output_dir=tmp_path / "default",
        max_records=3,
        library_dir=tmp_path / "empty",
        progress=False,
    )
    manifest = json.loads(result["manifest_path"].read_text())
    assert manifest["build_execution"] == "parallel-shards-v1"
    assert manifest["phase_metrics"]["prepare_seconds"] > 0
    assert manifest["payload_origins"] == {"items": "parsed_record"}
    assert result["resolution_stats"]["scope"] == "unique_entity_keys"
    assert not multiprocessing.active_children()


def test_zero_workers_cannot_select_another_pipeline(tmp_path, dataset):
    with pytest.raises(ValueError, match="At least one preparation worker"):
        build_resource("fixture", version="1", output_dir=tmp_path / "zero", batch_workers=0)
    assert not (tmp_path / "zero/resources/fixture/1").exists()


def test_preparation_processes_are_reused_across_batches(tmp_path, dataset, monkeypatch):
    from multiprocessing.process import BaseProcess

    original = BaseProcess.start
    started = []

    def start(process):
        original(process)
        started.append(process.pid)

    monkeypatch.setattr(BaseProcess, "start", start)
    result = build_resource(
        "fixture",
        version="1",
        output_dir=tmp_path / "reuse",
        batch_workers=2,
        max_records=12,
        max_batch_records=1,
        library_dir=tmp_path / "empty",
        progress=False,
    )
    assert result["batch_metrics"]["flushes"] == 11  # one intentionally skipped mapper result
    assert len(started) == 2
    assert not multiprocessing.active_children()


def test_mapper_ignores_module_attribute_named_annotations():
    from omnipath_build.two_phase import load_mapper
    from pypath.inputs_v2 import hpo

    assert load_mapper("pypath.inputs_v2.hpo", "annotations") is hpo.annotations_schema
