"""Replay every published source record through the native resource pipeline.

No downloads or reference builds. Original dataset:index IDs and exact stored JSON
text survive replay; rows discarded before publication cannot be reconstructed.
Run serially in a new private output root, then review artifacts before promotion.
"""

from __future__ import annotations

import argparse
from collections import Counter
from contextlib import contextmanager, ExitStack
import hashlib
import json
import os
from pathlib import Path
import pickle
import socket
import sys
import tempfile
from time import perf_counter
from unittest.mock import patch

import duckdb
import pyarrow.parquet as pq


ROW_ID_PATTERN = r"^([A-Za-z0-9_.-]+):([0-9]+)$"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_input(directory):
    """Attest the source payload file and recover the original dataset scope."""
    from omnipath_core.versioning import validate_source, validate_version
    from omnipath_core.schema import PAYLOAD_SCHEMA

    directory = Path(directory).resolve(strict=True)
    manifest_path = directory / "build_manifest.json"
    manifest_bytes = manifest_path.read_bytes()
    manifest = json.loads(manifest_bytes)
    source = validate_source(manifest["resource"])
    version = validate_version(manifest["version"])
    explicit = manifest.get("datasets")
    origins = manifest.get("payload_origins")
    require(
        explicit is None or isinstance(explicit, list), "Manifest datasets must be a list or null"
    )
    require(
        origins is None or isinstance(origins, dict), "Manifest payload_origins must be an object"
    )
    datasets = explicit if explicit is not None else list(origins or {})
    require(
        all(isinstance(d, str) and d for d in datasets) and len(set(datasets)) == len(datasets),
        "Original manifest dataset scope must contain unique nonempty names",
    )
    if origins:
        require(set(origins) == set(datasets), "Manifest datasets and payload_origins disagree")
        require(
            all(value == "parsed_record" for value in origins.values()),
            "Mapped-record fallback payloads cannot be replayed through raw input mappers",
        )
    payloads = directory / "evidence_payloads.parquet"
    expected = manifest["files"][payloads.name]
    before = payloads.stat()
    require(pq.read_schema(payloads).equals(PAYLOAD_SCHEMA), "Unsupported published payload schema")
    require(pq.read_metadata(payloads).num_rows == expected["rows"], "Payload row count changed")
    require(before.st_size == expected["size_bytes"], "Payload file size changed")
    checksum = digest_file(payloads)
    after = payloads.stat()
    require(
        (before.st_size, before.st_mtime_ns, before.st_ctime_ns)
        == (after.st_size, after.st_mtime_ns, after.st_ctime_ns),
        "Published input changed during verification",
    )
    require(checksum == expected["sha256"], "Payload checksum differs from original manifest")
    return {
        "source": source,
        "version": version,
        "datasets": datasets,
        "dataset_scope_origin": "manifest_datasets"
        if explicit
        else "manifest_payload_origins"
        if origins
        else "published_row_ids_plus_current_discovery",
        "payloads_path": payloads,
        "payloads_sha256": checksum,
        "manifest_path": manifest_path,
        "manifest_sha256": hashlib.sha256(manifest_bytes).hexdigest(),
        "reference": manifest.get("reference"),
    }


class ReplayRows:
    """Disk-backed deduplication and bounded source iteration."""

    def __init__(self, payloads, source, datasets, workdir, *, memory_limit="256MB", cap=None):
        require(cap is None or cap > 0, "Optional record cap must be positive")
        self.source, self.datasets, self.cap = source, list(datasets), cap
        self.consumed = Counter({dataset: 0 for dataset in datasets})
        self.current = None
        self.workdir = Path(workdir)
        self.started_datasets = set()
        self.conn = duckdb.connect(str(self.workdir / "replay.duckdb"))
        try:
            self.conn.execute("SET memory_limit=?", [memory_limit])
            self.conn.execute("SET threads=1")
            self.conn.execute("SET preserve_insertion_order=false")
            self.conn.execute("SET temp_directory=?", [str(Path(workdir) / "spill")])
            bad = self.conn.execute(
                """SELECT source,row_id FROM read_parquet(?)
                   WHERE source IS DISTINCT FROM ? OR row_id IS NULL OR payload_json IS NULL
                   OR NOT regexp_full_match(row_id, ?) LIMIT 1""",
                [str(payloads), source, ROW_ID_PATTERN],
            ).fetchone()
            require(bad is None, f"Unscoped, invalid or wrong-source published record: {bad}")
            # Exact strings, rather than semantic JSON equality, determine conflicts.
            self.conn.execute(
                """CREATE TABLE replay_rows AS
                   SELECT row_id, regexp_extract(row_id, ?, 1) AS dataset,
                          regexp_extract(row_id, ?, 2) AS original_index,
                          min(payload_json) AS payload_json,
                          min(payload_json) <> max(payload_json) AS conflict,
                          count(*) AS copies
                   FROM read_parquet(?) GROUP BY row_id""",
                [ROW_ID_PATTERN, ROW_ID_PATTERN, str(payloads)],
            )
            conflict = self.conn.execute(
                "SELECT row_id FROM replay_rows WHERE conflict LIMIT 1"
            ).fetchone()
            require(
                conflict is None, f"Conflicting exact raw variants for original record: {conflict}"
            )
            unknown = self.conn.execute(
                "SELECT DISTINCT dataset FROM replay_rows WHERE dataset NOT IN (SELECT unnest(?))",
                [self.datasets],
            ).fetchall()
            require(not unknown, f"Payload dataset absent from original manifest scope: {unknown}")
            self.counts = dict.fromkeys(datasets, 0)
            self.counts.update(
                dict(
                    self.conn.execute(
                        "SELECT dataset,count(*) FROM replay_rows GROUP BY dataset"
                    ).fetchall()
                )
            )
            self.selected_counts = {
                dataset: min(count, cap) if cap is not None else count
                for dataset, count in self.counts.items()
            }
            self.payload_rows = self.conn.execute(
                "SELECT coalesce(sum(copies),0) FROM replay_rows"
            ).fetchone()[0]
        except BaseException:
            self.close()
            raise

    def iterator(self, dataset):
        def raw(**_kwargs):
            require(dataset not in self.started_datasets, f"Dataset replayed twice: {dataset}")
            self.started_datasets.add(dataset)
            name = hashlib.sha256(dataset.encode()).hexdigest()[:24]
            spool = self.workdir / f"ordered-{name}.parquet"
            try:
                query = """SELECT original_index,payload_json FROM replay_rows WHERE dataset=$dataset
                    ORDER BY length(ltrim(original_index,'0')),ltrim(original_index,'0'),original_index"""
                parameters = {"dataset": dataset, "spool": str(spool)}
                if self.cap is not None:
                    query += " LIMIT $cap"
                    parameters["cap"] = self.cap
                # A COPY sink avoids materializing the full raw SQL result in
                # the connection before fetchmany. The sort can spill privately.
                self.conn.execute(
                    "COPY (" + query + ") TO $spool (FORMAT PARQUET, COMPRESSION ZSTD, "
                    "ROW_GROUP_SIZE_BYTES '8MB')",
                    parameters,
                )
                parquet = pq.ParquetFile(
                    spool, memory_map=False, pre_buffer=False, buffer_size=64 * 1024
                )
                try:
                    batches = parquet.iter_batches(batch_size=32, use_threads=False)
                    try:
                        for batch in batches:
                            for index in range(batch.num_rows):
                                item = batch.slice(index, 1).to_pylist()[0]
                                original_index, payload = (
                                    item["original_index"],
                                    item["payload_json"],
                                )
                                record = json.loads(payload)
                                require(
                                    isinstance(record, dict),
                                    f"Raw source object required: {dataset}:{original_index}",
                                )
                                self.current = (dataset, original_index, record, payload)
                                self.consumed[dataset] += 1
                                yield record
                            del batch
                    finally:
                        batches.close()
                finally:
                    parquet.close()
            finally:
                self.current = None
                spool.unlink(missing_ok=True)

        return raw

    def selected_sql(self):
        if self.cap is None:
            return "SELECT row_id,payload_json,dataset FROM replay_rows"
        return """SELECT row_id,payload_json,dataset FROM replay_rows
            QUALIFY row_number() OVER (PARTITION BY dataset
                ORDER BY length(ltrim(original_index,'0')),ltrim(original_index,'0'),original_index) <= ?"""

    def verify_output(self, result):
        require(
            dict(self.consumed) == self.selected_counts,
            "Native build did not consume selected rows",
        )
        parameters = [] if self.cap is None else [self.cap]
        self.conn.execute("CREATE TEMP TABLE selected_rows AS " + self.selected_sql(), parameters)
        self.conn.execute(
            """CREATE TEMP TABLE output_rows AS
               SELECT row_id,source,min(payload_json) AS payload_json,
                      min(payload_json) <> max(payload_json) AS conflict,
                      bool_or(payload_json IS NULL) AS missing_payload
               FROM read_parquet(?) GROUP BY row_id,source""",
            [str(result["files"]["evidence_payloads"])],
        )
        altered = self.conn.execute(
            """SELECT o.row_id FROM output_rows o LEFT JOIN selected_rows s USING (row_id)
               WHERE s.row_id IS NULL OR o.source IS DISTINCT FROM ? OR o.conflict
                   OR o.missing_payload OR o.payload_json IS DISTINCT FROM s.payload_json LIMIT 1""",
            [self.source],
        ).fetchone()
        require(altered is None, f"Replayed record identity or exact raw text changed: {altered}")
        missing = self.conn.execute(
            """SELECT s.row_id FROM selected_rows s LEFT JOIN output_rows o USING (row_id)
               WHERE o.row_id IS NULL LIMIT 1"""
        ).fetchone()
        require(
            missing is None,
            f"Current mapper discarded a previously published source record: {missing}",
        )
        return {
            "published_payload_rows": self.payload_rows,
            "published_source_records": sum(self.counts.values()),
            "replayed_source_records": sum(self.consumed.values()),
            "datasets": {
                d: {"published_records": self.counts[d], "replayed_records": self.consumed[d]}
                for d in self.datasets
            },
            "exact_original_ids_and_payloads_verified": True,
        }

    def write_source_index(self, path):
        """Durable per-record hashes without retaining records in Python memory."""
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        self.conn.execute(
            """COPY (SELECT dataset,row_id,sha256(payload_json) AS source_record_sha256
               FROM selected_rows) TO ? (FORMAT PARQUET, COMPRESSION ZSTD)""",
            [str(path)],
        )
        return {"path": str(path), "rows": sum(self.consumed.values()), "sha256": digest_file(path)}

    def close(self):
        self.conn.close()


class ParentInputPickle:
    """Adapt only the native parent's input tuple; workers receive ordinary data."""

    def __init__(self, rows):
        self.rows = rows

    def dump(self, value, handle, *args, **kwargs):
        require(
            isinstance(value, tuple)
            and len(value) == 3
            and isinstance(value[0], int)
            and self.rows.current is not None,
            "Unexpected native parent input tuple; review replay adapter compatibility",
        )
        _dataset, original_index, record, original_payload = self.rows.current
        require(value[1] is record, "Native input does not match the active original source record")
        return pickle.dump((original_index, record, original_payload), handle, *args, **kwargs)

    def __getattr__(self, name):
        return getattr(pickle, name)


def blocked_download(*_args, **_kwargs):
    raise RuntimeError("Source downloads are disabled during published-resource replay")


class OfflineCurlHandle:
    """Delegate local curl configuration, rejecting every transfer method."""

    def __init__(self, handle):
        self.handle = handle

    def __getattr__(self, name):
        if name.startswith("perform"):
            return blocked_download
        return getattr(self.handle, name)

    def duphandle(self):
        return OfflineCurlHandle(self.handle.duphandle())

    def __enter__(self):
        return self

    def __exit__(self, *_args):
        self.handle.close()


@contextmanager
def offline_network():
    """Block Python sockets and known native download transports in every process."""
    from dlmachine import DownloadManager, _downloader

    connect, connect_ex = socket.socket.connect, socket.socket.connect_ex

    def guarded(method):
        def call(sock, *args, **kwargs):
            if sock.family in (socket.AF_INET, socket.AF_INET6):
                return blocked_download()
            return method(sock, *args, **kwargs)

        return call

    with ExitStack() as stack:
        stack.enter_context(patch.object(socket.socket, "connect", guarded(connect)))
        stack.enter_context(patch.object(socket.socket, "connect_ex", guarded(connect_ex)))
        stack.enter_context(patch.object(DownloadManager, "download", blocked_download))
        for backend in (_downloader.CurlDownloader, _downloader.RequestsDownloader):
            # Covers instances created before the guard, including native handles.
            stack.enter_context(patch.object(backend, "download", blocked_download))
        # Existing imported aliases still reach the patched DownloadManager or
        # backend. Patch loaded helper bindings as well to reject calls earlier.
        for name in ("pypath.share.downloads", "pypath.inputs_v2.base"):
            module = sys.modules.get(name)
            if module is not None:
                stack.enter_context(patch.object(module, "download_and_open", blocked_download))
        legacy = sys.modules.get("pypath.share.curl")
        if legacy is not None:
            stack.enter_context(patch.object(legacy.Curl, "curl_call", blocked_download))
        try:
            import pycurl
        except ModuleNotFoundError as exc:
            if exc.name != "pycurl":
                raise
        else:
            native_curl = pycurl.Curl

            def local_handle(*args, **kwargs):
                return OfflineCurlHandle(native_curl(*args, **kwargs))

            stack.enter_context(patch.object(pycurl, "Curl", local_handle))
            if hasattr(pycurl, "CurlMulti"):
                # Native multi transfers can bypass Curl.perform entirely.
                stack.enter_context(patch.object(pycurl, "CurlMulti", blocked_download))
        yield


def offline_prepare_worker(*args, **kwargs):
    from omnipath_build.two_phase import prepare_worker

    with offline_network():
        prepare_worker(*args, **kwargs)


def replay_resource(
    original,
    root,
    version,
    library,
    *,
    batch_workers=1,
    resource_ram_bytes=2 * 1024**3,
    memory_limit="256MB",
    max_records=None,
    min_free_disk_bytes=1024**3,
):
    from omnipath_build.discovery import discover_datasets
    from omnipath_build.pipeline import build_resource
    from omnipath_build import two_phase

    require(batch_workers > 0, "Preparation workers must be positive")
    started = perf_counter()
    root = Path(root)
    cache = root / "pypath-data"
    with offline_network():
        _, datasets, _ = discover_datasets(
            original["source"], datasets=original["datasets"] or None, cache_dir=cache
        )
    dataset_names = original["datasets"] or [dataset.dataset_name for dataset in datasets]
    with tempfile.TemporaryDirectory(prefix=".replay-", dir=root) as work:
        rows = ReplayRows(
            original["payloads_path"],
            original["source"],
            dataset_names,
            work,
            memory_limit=memory_limit,
            cap=max_records,
        )
        try:
            with offline_network(), ExitStack() as stack:
                for dataset in datasets:
                    require(
                        dataset.raw_dataset is not None and dataset.mapper is not None,
                        f"Native raw mapper required for {dataset.dataset_name}",
                    )
                    stack.enter_context(
                        patch.object(
                            dataset.raw_dataset, "raw", rows.iterator(dataset.dataset_name)
                        )
                    )
                    # A replay reproduces what was published, an empty dataset included.
                    stack.enter_context(patch.object(dataset.raw_dataset, "allow_empty", True))
                stack.enter_context(patch.object(two_phase, "pickle", ParentInputPickle(rows)))
                stack.enter_context(
                    patch.object(two_phase, "prepare_worker", offline_prepare_worker)
                )
                result = build_resource(
                    original["source"],
                    version=version,
                    output_dir=root,
                    library_dir=library,
                    cache_dir=cache,
                    datasets=dataset_names,
                    max_records=max_records,
                    batch_workers=batch_workers,
                    resource_ram_bytes=resource_ram_bytes,
                    min_free_disk_bytes=min_free_disk_bytes,
                    max_batch_records=4096,
                    batch_size=2000,
                    max_batch_relations=2000,
                    max_batch_bytes=32 * 1024**2,
                    raw_batch_bytes=16 * 1024**2,
                    progress=True,
                )
            coverage = rows.verify_output(result)
            invalid_record = rows.conn.execute(
                """SELECT r.relation_key,ev.row_id FROM read_parquet(?) r
                   JOIN read_parquet(?) ev USING (relation_id)
                   LEFT JOIN selected_rows s ON s.row_id=ev.row_id
                   WHERE r.statement_kind='relation' AND r.subject_type='molecular_activity'
                     AND r.predicate IN ('has_input','has_output','enabled_by')
                     AND (s.row_id IS NULL OR ev.source IS DISTINCT FROM ?
                       OR ev.dataset IS DISTINCT FROM s.dataset)
                   LIMIT 1""",
                [
                    str(result["files"]["relation"]),
                    str(result["files"]["relation_evidence"]),
                    original["source"],
                ],
            ).fetchone()
            require(
                invalid_record is None,
                f"Reaction occurrence lost its original source record: {invalid_record}",
            )
            coverage["source_record_index"] = rows.write_source_index(
                root / "replay_source_records" / f"{original['source']}.parquet"
            )
            require(
                digest_file(original["payloads_path"]) == original["payloads_sha256"],
                "Original published payload file changed during replay",
            )
        finally:
            rows.close()
    require(
        digest_file(original["manifest_path"]) == original["manifest_sha256"],
        "Original manifest changed during replay",
    )
    return {
        "resource": original["source"],
        "original_manifest_path": str(original["manifest_path"]),
        "original_payloads_path": str(original["payloads_path"]),
        "dataset_scope_origin": original["dataset_scope_origin"],
        "published_datasets": sorted(
            dataset for dataset, item in coverage["datasets"].items() if item["published_records"]
        ),
        "original_version": original["version"],
        "version": version,
        "original_manifest_sha256": original["manifest_sha256"],
        "original_payloads_sha256": original["payloads_sha256"],
        "original_reference": original["reference"],
        "manifest_path": str(result["manifest_path"]),
        "seconds": round(perf_counter() - started, 6),
        "phase_metrics": result["phase_metrics"],
        **coverage,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-resource-dir", action="append", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--version", required=True)
    parser.add_argument("--library-dir", required=True, type=Path)
    parser.add_argument("--human-gem-xrefs", type=Path)
    parser.add_argument("--batch-workers", type=int, default=1)
    parser.add_argument("--resource-ram-gb", type=int, default=2)
    parser.add_argument("--duckdb-threads", type=int, default=1)
    parser.add_argument(
        "--memory-limit", default="256MB", help="Replay deduplication memory, with disk spill"
    )
    parser.add_argument(
        "--max-records", type=int, help="Optional per-dataset test cap; absent means all rows"
    )
    parser.add_argument("--min-free-disk-gb", type=int, default=1)
    args = parser.parse_args(argv)
    from omnipath_resolver.canonical.library import pin_library
    from omnipath_core.versioning import validate_version

    version = validate_version(args.version)
    originals = [read_input(path) for path in args.input_resource_dir]
    require(
        len({item["source"] for item in originals}) == len(originals), "Duplicate resource input"
    )
    require(
        args.batch_workers > 0
        and args.resource_ram_gb > 0
        and args.duckdb_threads > 0
        and args.min_free_disk_gb >= 0,
        "Worker, memory, thread and disk budgets are invalid",
    )
    library = pin_library(args.library_dir.resolve(strict=True))
    require(
        (library / "manifest.json").is_file(),
        "An existing immutable reference generation is required",
    )
    root = args.output_dir.resolve()
    require(not root.exists(), "Use a new private output directory")
    require(not library.is_relative_to(root), "Output must not contain the reference")
    require(
        all(not item["payloads_path"].is_relative_to(root) for item in originals),
        "Output must not contain original inputs",
    )
    if any(item["source"] == "metatlas" for item in originals):
        xrefs = args.human_gem_xrefs or os.environ.get("OMNIPATH_HUMAN_GEM_XREFS")
        require(xrefs, "Metatlas replay requires an existing --human-gem-xrefs TSV; no downloads")
        os.environ["OMNIPATH_HUMAN_GEM_XREFS"] = str(Path(xrefs).resolve(strict=True))
    root.mkdir(parents=True)
    reference_checksum = digest_file(library / "manifest.json")
    report = {
        "status": "running",
        "version": version,
        "reference": str(library),
        "reference_manifest_sha256": reference_checksum,
        "max_records_per_dataset": args.max_records,
        "coverage": "every published original record; unpublished discarded source rows unavailable",
        "replay_script_sha256": digest_file(Path(__file__)),
        "preparation_workers": args.batch_workers,
        "resource_ram_bytes": args.resource_ram_gb * 1024**3,
        "deduplication_memory_limit": args.memory_limit,
        "duckdb_threads": args.duckdb_threads,
        "resources": {},
    }
    try:
        with patch.dict(os.environ, {"OMNIPATH_BUILD_DUCKDB_THREADS": str(args.duckdb_threads)}):
            for original in originals:
                resource = replay_resource(
                    original,
                    root,
                    version,
                    library,
                    batch_workers=args.batch_workers,
                    resource_ram_bytes=args.resource_ram_gb * 1024**3,
                    memory_limit=args.memory_limit,
                    max_records=args.max_records,
                    min_free_disk_bytes=args.min_free_disk_gb * 1024**3,
                )
                report["resources"][original["source"]] = resource
                require(
                    digest_file(library / "manifest.json") == reference_checksum,
                    "Reference changed during replay",
                )
                print(json.dumps(resource), flush=True)
                (root / "replay_report.json").write_text(json.dumps(report, indent=2) + "\n")
        report["status"] = "success"
    except BaseException as exc:
        report.update(status="failed", error_type=type(exc).__name__, error=str(exc))
        raise
    finally:
        (root / "replay_report.json").write_text(json.dumps(report, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
