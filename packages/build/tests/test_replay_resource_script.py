"""Published source replay preserves exact occurrence identity and stays offline."""

import hashlib
import io
import json
import pickle
import socket
import sys
from types import ModuleType, SimpleNamespace

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from omnipath_core.schema import PAYLOAD_SCHEMA
from scripts import replay_resource as replay


def write_payloads(path, records):
    pq.write_table(pa.Table.from_pylist(records, schema=PAYLOAD_SCHEMA), path)
    return path


def payload(row_id="reactions:7", raw='{ "value": "μ" }', source="rhea"):
    return {
        "source": source,
        "row_id": row_id,
        "payload_json": raw,
        "entity_key": "E1",
        "relation_key": None,
    }


def published(tmp_path, records, *, datasets=None, origins=None, source="rhea"):
    directory = tmp_path / "original"
    directory.mkdir()
    path = write_payloads(directory / "evidence_payloads.parquet", records)
    manifest = {
        "resource": source,
        "version": "2026.09.1",
        "datasets": datasets,
        "payload_origins": origins,
        "files": {
            path.name: {
                "rows": len(records),
                "size_bytes": path.stat().st_size,
                "sha256": replay.digest_file(path),
            }
        },
    }
    (directory / "build_manifest.json").write_text(json.dumps(manifest))
    return directory, manifest


def rows_at(tmp_path, records, datasets=("reactions",), cap=None):
    path = write_payloads(tmp_path / "raw.parquet", records)
    return replay.ReplayRows(path, "rhea", datasets, tmp_path, cap=cap)


def test_manifest_retains_zero_output_datasets_and_attests_payload(tmp_path):
    directory, manifest = published(
        tmp_path,
        [payload()],
        origins={"reactions": "parsed_record", "transport_reactions": "parsed_record"},
    )
    original = replay.read_input(directory)
    assert original["datasets"] == ["reactions", "transport_reactions"]
    assert original["payloads_sha256"] == manifest["files"]["evidence_payloads.parquet"]["sha256"]
    assert original["manifest_sha256"] == replay.digest_file(directory / "build_manifest.json")


@pytest.mark.parametrize(
    "datasets,origins,message",
    [
        (["reactions"], {"other": "parsed_record"}, "disagree"),
        (["reactions"], {"reactions": "mapped_record_fallback"}, "fallback"),
        (["reactions", "reactions"], None, "unique nonempty"),
    ],
)
def test_manifest_rejects_ambiguous_or_unreplayable_scope(tmp_path, datasets, origins, message):
    directory, _ = published(tmp_path, [payload()], datasets=datasets, origins=origins)
    with pytest.raises(ValueError, match=message):
        replay.read_input(directory)


def test_manifest_rejects_modified_input(tmp_path):
    directory, manifest = published(tmp_path, [payload()], datasets=["reactions"])
    manifest["files"]["evidence_payloads.parquet"]["sha256"] = "0" * 64
    (directory / "build_manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="checksum"):
        replay.read_input(directory)


def test_deduplicates_owners_and_preserves_gaps_exact_text_and_numeric_order(tmp_path):
    records = [payload("reactions:10"), payload("reactions:2"), payload("reactions:2")]
    records[2]["entity_key"] = None
    records[2]["relation_key"] = "R1"
    rows = rows_at(tmp_path, records, datasets=("reactions", "transport_reactions"))
    try:
        raw = rows.iterator("reactions")()
        adapter = replay.ParentInputPickle(rows)
        packets = io.BytesIO()
        for index, record in enumerate(raw):
            adapter.dump((index, record, json.dumps(record)), packets)
        packets.seek(0)
        assert pickle.load(packets) == ("2", {"value": "μ"}, records[0]["payload_json"])
        assert pickle.load(packets)[0] == "10"
        assert rows.counts == {"reactions": 2, "transport_reactions": 0}
        assert rows.payload_rows == 3
        assert list(rows.iterator("transport_reactions")()) == []
        assert rows.selected_counts == {"reactions": 2, "transport_reactions": 0}
    finally:
        rows.close()


def test_rejects_conflicting_exact_variants_even_for_equivalent_json(tmp_path):
    with pytest.raises(ValueError, match="Conflicting exact raw"):
        rows_at(tmp_path, [payload(raw='{"value":1}'), payload(raw='{ "value": 1 }')])


@pytest.mark.parametrize(
    "record",
    [
        payload(source="other"),
        payload(source=None),
        payload(row_id=None),
        payload(row_id="7"),
        payload(row_id="reactions:-1"),
        payload(row_id="reactions:x"),
        payload(row_id="other:7"),
        payload(raw=None),
    ],
)
def test_rejects_wrong_source_invalid_ids_and_unknown_datasets(tmp_path, record):
    with pytest.raises(ValueError):
        rows_at(tmp_path, [record])


def test_cap_is_optional_and_only_limits_selected_unique_source_rows(tmp_path):
    rows = rows_at(tmp_path, [payload(f"reactions:{i}") for i in (1, 9, 12)], cap=2)
    try:
        assert len(list(rows.iterator("reactions")())) == 2
        assert rows.counts["reactions"] == 3
        assert rows.consumed["reactions"] == 2
    finally:
        rows.close()


def test_raw_nonobject_is_rejected_without_network_or_mapper_changes(tmp_path):
    rows = rows_at(tmp_path, [payload(raw="[1,2]")])
    try:
        with pytest.raises(ValueError, match="Raw source object"):
            list(rows.iterator("reactions")())
    finally:
        rows.close()


def test_adapter_changes_only_expected_active_input_tuple(tmp_path):
    rows = rows_at(tmp_path, [payload()])
    try:
        raw = rows.iterator("reactions")()
        record = next(raw)
        adapter = replay.ParentInputPickle(rows)
        for packet in ((0, {}, "{}"), ("0", record, "{}"), (0, record), None):
            with pytest.raises(ValueError):
                adapter.dump(packet, io.BytesIO())
        assert adapter.HIGHEST_PROTOCOL == pickle.HIGHEST_PROTOCOL
        raw.close()
        with pytest.raises(ValueError):
            adapter.dump((0, record, "{}"), io.BytesIO())
    finally:
        rows.close()


def test_socket_guard_denies_internet_and_restores_on_exception():
    original = socket.socket.connect
    with pytest.raises(RuntimeError, match="downloads are disabled"):
        with replay.offline_network():
            with socket.socket(socket.AF_INET) as sock:
                sock.connect(("127.0.0.1", 9))
    assert socket.socket.connect is original


def test_real_native_mapper_worker_preserves_original_ids_hashes_and_zero_dataset(
    tmp_path, monkeypatch
):
    from omnipath_build import two_phase
    from pypath.inputs_v2 import tcdb

    raw = {
        "transporter_uniprot": "P11166",
        "substrate_id": "CHEBI:4167",
        "substrate_name": "glucose",
        "tcid": "2.A.1.1.28",
    }
    texts = [json.dumps(raw, ensure_ascii=False, indent=2), json.dumps({**raw, "note": "μ"})]
    records = [payload(f"transports:{i}", text, source="tcdb") for i, text in zip((2, 9), texts)]
    directory, _ = published(
        tmp_path,
        records,
        origins={"transporters": "parsed_record", "transports": "parsed_record"},
        source="tcdb",
    )
    original = replay.read_input(directory)
    output = tmp_path / "output"
    output.mkdir()
    native_raw = tcdb.resource.datasets()["transports"].raw
    native_mapper = tcdb.resource.datasets()["transports"].mapper
    native_pickle, native_worker = two_phase.pickle, two_phase.prepare_worker
    monkeypatch.setenv("OMNIPATH_BUILD_DUCKDB_THREADS", "1")
    report = replay.replay_resource(
        original, output, "2026.09.2", None, resource_ram_bytes=1024**3, min_free_disk_bytes=0
    )
    assert report["replayed_source_records"] == 2
    assert report["datasets"]["transporters"]["replayed_records"] == 0
    assert tcdb.resource.datasets()["transports"].raw == native_raw
    assert tcdb.resource.datasets()["transports"].mapper is native_mapper
    assert not tcdb.resource.datasets()["transporters"].allow_empty
    assert two_phase.pickle is native_pickle
    assert two_phase.prepare_worker is native_worker
    compiled = output / "resources" / "tcdb" / "2026.09.2"
    evidence = pq.read_table(compiled / "relation_evidence.parquet").to_pylist()
    assert {ev["row_id"] for ev in evidence} == {"transports:2", "transports:9"}
    expected = {
        f"transports:{i}": hashlib.sha256(text.encode()).hexdigest()
        for i, text in zip((2, 9), texts)
    }
    index = pq.read_table(report["source_record_index"]["path"]).to_pylist()
    assert {(row["row_id"], row["source_record_sha256"]) for row in index} == set(expected.items())
    manifest = json.loads((compiled / "build_manifest.json").read_text())
    assert sorted(manifest["datasets"]) == ["transporters", "transports"]
    assert manifest["max_records"] is None


def test_native_parent_bindings_restore_when_build_raises(tmp_path, monkeypatch):
    from omnipath_build import pipeline, two_phase
    from pypath.inputs_v2 import rhea

    directory, _ = published(tmp_path, [payload()], datasets=["reactions"])
    original = replay.read_input(directory)
    output = tmp_path / "output"
    output.mkdir()
    native_pickle, native_worker = two_phase.pickle, two_phase.prepare_worker
    native_raw = rhea.resource.datasets()["reactions"].raw

    def fail(*_args, **_kwargs):
        assert two_phase.pickle is not native_pickle
        assert two_phase.prepare_worker is replay.offline_prepare_worker
        raise RuntimeError("deliberate native failure")

    monkeypatch.setattr(pipeline, "build_resource", fail)
    with pytest.raises(RuntimeError, match="deliberate native failure"):
        replay.replay_resource(original, output, "2", None)
    assert two_phase.pickle is native_pickle
    assert two_phase.prepare_worker is native_worker
    assert rhea.resource.datasets()["reactions"].raw == native_raw
    assert not list(output.glob(".replay-*"))


def test_legacy_empty_scope_replays_published_rows_and_infers_current_zero_datasets(
    tmp_path, monkeypatch
):
    from pypath.inputs_v2 import tcdb

    raw = {"transporter_uniprot": "P11166", "substrate_id": "CHEBI:4167", "tcid": "2.A.1.1.28"}
    directory, _ = published(
        tmp_path, [payload("transports:8", json.dumps(raw), source="tcdb")], origins={}, source="tcdb"
    )
    original = replay.read_input(directory)
    assert original["datasets"] == []
    assert original["dataset_scope_origin"] == "published_row_ids_plus_current_discovery"
    output = tmp_path / "output"
    output.mkdir()
    monkeypatch.setenv("OMNIPATH_BUILD_DUCKDB_THREADS", "1")
    report = replay.replay_resource(
        original, output, "2", None, resource_ram_bytes=1024**3, min_free_disk_bytes=0
    )
    assert report["replayed_source_records"] == 1
    assert report["published_datasets"] == ["transports"]
    assert report["dataset_scope_origin"] == "published_row_ids_plus_current_discovery"
    assert set(report["datasets"]) == set(tcdb.resource.datasets())
    assert report["datasets"]["transporters"]["replayed_records"] == 0


def test_mixed_valid_and_null_output_payloads_are_rejected(tmp_path):
    records = [payload()]
    rows = rows_at(tmp_path, records)
    try:
        list(rows.iterator("reactions")())
        null_copy = {
            **records[0],
            "entity_key": None,
            "relation_key": "other owner",
            "payload_json": None,
        }
        output = write_payloads(tmp_path / "output.parquet", [records[0], null_copy])
        with pytest.raises(ValueError, match="exact raw text changed"):
            rows.verify_output({"files": {"evidence_payloads": output}})
    finally:
        rows.close()


def test_native_curl_transport_proxy_rejects_perform_variants_and_duplicate_handles(monkeypatch):
    fake = ModuleType("pycurl")
    called = []

    class NativeHandle:
        def perform(self):
            called.append("native perform")

        perform_rb = perform_rs = perform

        def duphandle(self):
            return NativeHandle()

        def close(self):
            called.append("close")

        def setopt(self, *_args):
            called.append("local configuration")

    fake.Curl = NativeHandle
    fake.CurlMulti = lambda: called.append("native multi")
    monkeypatch.setitem(sys.modules, "pycurl", fake)
    original = fake.Curl
    with replay.offline_network():
        with fake.Curl() as curl:
            curl.setopt("local option", "value")
            for method in (
                curl.perform,
                curl.perform_rb,
                curl.perform_rs,
                curl.duphandle().perform,
                fake.CurlMulti,
            ):
                with pytest.raises(RuntimeError, match="downloads are disabled"):
                    method()
    assert called == ["local configuration", "close"]
    assert fake.Curl is original


def test_real_download_backend_and_preexisting_native_handle_entrypoints_are_blocked(monkeypatch):
    from dlmachine import DownloadManager, _downloader

    legacy = ModuleType("pypath.share.curl")

    class LegacyCurl:
        def curl_call(self):
            pytest.fail("Legacy transport entrypoint executed")

    legacy.Curl = LegacyCurl
    monkeypatch.setitem(sys.modules, "pypath.share.curl", legacy)
    download = _downloader.CurlDownloader.download
    called = []
    preexisting = object.__new__(_downloader.CurlDownloader)
    preexisting._show_progress = False
    preexisting._progress_bar = None
    preexisting.handler = SimpleNamespace(perform=lambda: called.append("native perform"))
    with replay.offline_network():
        for call in (
            preexisting.download,
            LegacyCurl().curl_call,
            lambda: DownloadManager.download(None, "https://example.invalid"),
            lambda: _downloader.RequestsDownloader.download(None),
        ):
            with pytest.raises(RuntimeError, match="downloads are disabled"):
                call()
    assert called == []
    assert _downloader.CurlDownloader.download is download


def test_real_pycurl_transfer_is_rejected_without_invoking_native_transport():
    pycurl = pytest.importorskip("pycurl")
    original = pycurl.Curl
    with replay.offline_network():
        with pycurl.Curl() as curl:
            with pytest.raises(RuntimeError, match="downloads are disabled"):
                curl.perform()
    assert pycurl.Curl is original


def test_offline_guard_preserves_direct_local_file_reads(tmp_path):
    path = tmp_path / "local.tsv"
    path.write_text("key\tvalue\n")
    with replay.offline_network():
        assert path.read_text() == "key\tvalue\n"


def test_offline_guard_preserves_cached_human_gem_xref_mapper(tmp_path, monkeypatch):
    from pypath.inputs_v2 import metatlas
    from pypath.inputs_v2.parsers.metatlas import _metabolite_xrefs

    path = tmp_path / "Human-GEM-metabolites.tsv"
    path.write_text("metsNoComp\tmetSmiles\nM1\tO\n")
    monkeypatch.setenv("OMNIPATH_HUMAN_GEM_XREFS", str(path))
    _metabolite_xrefs.cache_clear()
    try:
        with replay.offline_network():
            mapped = metatlas._map_with_structures(
                metatlas.metabolites_schema, {"human_gem_metabolite_id": "M1"}
            )
        assert any(
            str(identifier.type) == "smiles" and identifier.value == "O"
            for identifier in mapped.identifiers
        )
    finally:
        _metabolite_xrefs.cache_clear()


def test_ordered_spool_has_buffered_reads_and_closes_on_early_stop(tmp_path, monkeypatch):
    rows = rows_at(tmp_path, [payload("reactions:10"), payload("reactions:2")])
    native_reader = pq.ParquetFile
    options, readers = [], []

    class Reader:
        def __init__(self, path, **kwargs):
            options.append(kwargs)
            self.reader = native_reader(path, **kwargs)
            self.closed = False
            readers.append(self)

        def iter_batches(self, **kwargs):
            assert kwargs == {"batch_size": 32, "use_threads": False}
            yield from self.reader.iter_batches(**kwargs)

        def close(self):
            self.closed = True
            self.reader.close()

    monkeypatch.setattr(replay.pq, "ParquetFile", Reader)
    iterator = rows.iterator("reactions")()
    try:
        assert next(iterator) == {"value": "μ"}
        assert rows.current[1] == "2"
        assert list(tmp_path.glob("ordered-*.parquet"))
        iterator.close()
        assert rows.current is None
        assert readers[0].closed
        assert not list(tmp_path.glob("ordered-*.parquet"))
        assert options == [{"memory_map": False, "pre_buffer": False, "buffer_size": 64 * 1024}]
    finally:
        iterator.close()
        rows.close()


def test_compressed_large_records_use_sql_sink_not_materialized_raw_result(tmp_path):
    import psutil

    process = psutil.Process()
    text = json.dumps({"value": "μ", "blob": "x" * (4 * 1024**2)}, ensure_ascii=False, indent=2)
    ids = (100, 2, 9, 20, 7, 3, 41, 8)
    path = tmp_path / "large.parquet"
    pq.write_table(
        pa.Table.from_pylist([payload(f"reactions:{i}", text) for i in ids], schema=PAYLOAD_SCHEMA),
        path,
        compression="zstd",
    )
    assert path.stat().st_size < 1024**2
    rows = replay.ReplayRows(path, "rhea", ["reactions"], tmp_path, memory_limit="128MB")
    original_connection = rows.conn
    statements = []

    class SinkOnly:
        def execute(self, query, parameters):
            assert query.startswith("COPY (")  # A raw SELECT.execute() regresses the memory fix.
            statements.append(query)
            return original_connection.execute(query, parameters)

        def close(self):
            original_connection.close()

    rows.conn = SinkOnly()
    iterator = rows.iterator("reactions")()
    try:
        before = process.memory_info().rss
        first = next(iterator)
        assert process.memory_info().rss - before < 192 * 1024**2
        assert first == json.loads(text)
        assert rows.current[1] == "2"
        assert rows.current[3] == text
        assert (
            hashlib.sha256(rows.current[3].encode()).hexdigest()
            == hashlib.sha256(text.encode()).hexdigest()
        )
        assert len(statements) == 1
    finally:
        iterator.close()
        rows.close()
    assert not list(tmp_path.glob("ordered-*.parquet"))
