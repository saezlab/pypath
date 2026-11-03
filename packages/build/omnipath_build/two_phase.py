"""Parallel private observation shards, followed by one finalizer.

The pipeline publishes a resource version only after these phases succeed. The reference library must stay
immutable for the run. Input pickle files are private locally generated artifacts.
"""

from pathlib import Path
import json
import os
import pickle
import time
import importlib
from functools import lru_cache
from .contracts import BuildResult, WorkerResult, TwoPhaseConfig
from .discovery import DiscoveredDataset
from .progress import ProgressCallback, CancelCheck
from collections.abc import Sequence
from .provenance import file_fingerprint


@lru_cache(maxsize=16)
def load_mapper(module: str, dataset: str):
    from pypath.inputs_v2.base import Resource

    mod = importlib.import_module(module)
    ds = getattr(mod, dataset, None)
    if not hasattr(getattr(ds, "_raw_dataset", ds), "mapper"):
        ds = None
        for value in vars(mod).values():
            if isinstance(value, Resource) and dataset in value.datasets():
                ds = value.datasets()[dataset]
                break
    if ds is None:
        raise ValueError(f"Cannot load mapper for {module}:{dataset}")
    raw = getattr(ds, "_raw_dataset", ds)
    return getattr(raw, "mapper", None)


def prepare_worker(worker_id, tasks, results, root, source, library_dir, memory_limit, limits=None):
    """One reusable resolver/writer per worker; no shared mutable database."""
    os.environ["OMNIPATH_BUILD_DUCKDB_THREADS"] = "1"
    from omnipath_resolver import EntityResolver
    from .silver import SilverExtractor
    from .writer import ParquetWriter

    limits = limits or {"entities": 5000, "relations": 5000, "bytes": 64 * 1024**2, "records": 8192}
    resolver = writer = None
    try:
        resolver = EntityResolver(library_dir, defer_aliases=True, memory_limit=memory_limit)
        writer = ParquetWriter(
            Path(root) / f"worker-{worker_id}", library_dir=library_dir, memory_limit=memory_limit
        )
        started = time.monotonic()
        batch_metrics = {
            "flushes": 0,
            "peak_records": 0,
            "peak_relations": 0,
            "peak_estimated_bytes": 0,
        }
        while True:
            task = tasks.get()
            if task is None:
                break
            path, module, dataset, mapped = task
            mapper = load_mapper(module, dataset) if mapped else None
            extractor = SilverExtractor(source, dataset)
            rows = 0

            def flush():
                nonlocal extractor
                if extractor.entities or extractor.relations:
                    batch_metrics["flushes"] += 1
                    for key, value in [
                        ("peak_records", extractor.record_count),
                        ("peak_relations", len(extractor.relations)),
                        ("peak_estimated_bytes", extractor.estimated_bytes),
                    ]:
                        batch_metrics[key] = max(batch_metrics[key], value)
                    writer.append_observations(extractor, resolver)
                    extractor = SilverExtractor(source, dataset)

            with open(path, "rb") as handle:
                while True:
                    try:
                        index, raw, payload = pickle.load(handle)
                    except EOFError:
                        break
                    try:
                        record = mapper(raw) if mapper else raw
                        if record is not None:
                            extractor.process_record(
                                record, raw, f"{dataset}:{index}", index, payload_json=payload
                            )
                    except Exception as exc:
                        # Mappers fail loudly; say which input row to inspect.
                        raise RuntimeError(
                            f"{source}:{dataset} input row {index} ({module}) failed: "
                            f"{type(exc).__name__}: {exc}"
                        ) from exc
                    rows += 1
                    if (
                        len(extractor.entities) >= limits["entities"]
                        or len(extractor.relations) >= limits["relations"]
                        or extractor.estimated_bytes >= limits["bytes"]
                        or extractor.record_count >= limits["records"]
                    ):
                        flush()
            flush()
            Path(path).unlink()
            results.put({"worker": worker_id, "batch": Path(path).name, "rows": rows})
        stats = resolver.resolution_stats()
        resolver.export_resolution_keys(writer.output_dir / "resolution_keys.parquet")
        resolver.close()
        resolver = None
        shard = writer.seal_observation_shard()
        writer = None
        results.put(
            {
                "worker": worker_id,
                "shard": shard,
                "resolution": stats,
                "seconds": time.monotonic() - started,
                "batch_metrics": batch_metrics,
            }
        )
    except BaseException:
        import traceback

        results.put({"worker": worker_id, "error": traceback.format_exc()})
        raise
    finally:
        if resolver:
            resolver.close()
        if writer:
            writer.abort()


def build_two_phase(
    source: str,
    datasets: Sequence[DiscoveredDataset],
    target: str | Path,
    *,
    config: TwoPhaseConfig,
    on_progress: ProgressCallback | None = None,
    should_cancel: CancelCheck | None = None,
) -> BuildResult:
    """Prepare private shards, reap all workers, then aggregate atomically upstream."""
    import multiprocessing as mp
    import queue
    import shutil
    import tempfile
    from itertools import islice
    from .progress import check_cancel
    from .silver import canonical_json
    from .writer import ParquetWriter
    from .pipeline import _notify

    workers, cpus, ram_bytes = config.workers, config.cpus, config.ram_bytes
    final_memory, library_dir = config.final_memory, config.library_dir
    final_threads = config.final_threads or cpus
    max_records, force_refresh = config.max_records, config.force_refresh
    max_batch_records, raw_batch_bytes = config.max_batch_records, config.raw_batch_bytes
    min_free_disk_bytes = config.min_free_disk_bytes
    entity_limit, relation_limit, byte_limit = (
        config.entity_limit,
        config.relation_limit,
        config.byte_limit,
    )
    started = time.monotonic()
    target = Path(target)
    target.mkdir(parents=True, exist_ok=True)

    def notify(stage, message, **counts):
        _notify(
            on_progress,
            prints=False,
            pipeline="resource",
            stage=stage,
            status="running",
            message=message,
            counts=counts,
        )

    context = mp.get_context("spawn")
    workers = max(1, min(workers, cpus, max(1, ram_bytes // (1024**3))))
    sql_memory = f"{max(64 * 1024**2, min(512 * 1024**2, ram_bytes // (workers * 8)))}B"
    limits = {
        "entities": entity_limit,
        "relations": relation_limit,
        "records": min(max_batch_records, 8192),
        "bytes": min(byte_limit, 64 * 1024**2, max(1024**2, ram_bytes // (workers * 16))),
    }
    total_rows = 0
    finished: dict[int, WorkerResult] = {}
    input_fingerprints = {}
    with tempfile.TemporaryDirectory(prefix="phase-", dir=target) as work:
        root = Path(work)
        tasks, results = context.Queue(maxsize=workers * 2), context.Queue()
        processes = [
            context.Process(
                target=prepare_worker,
                args=(
                    i,
                    tasks,
                    results,
                    work,
                    source,
                    str(library_dir) if library_dir else None,
                    sql_memory,
                    limits,
                ),
            )
            for i in range(workers)
        ]

        def drain():
            nonlocal total_rows
            check_cancel(should_cancel)
            if shutil.disk_usage(target).free < min_free_disk_bytes:
                raise RuntimeError("Free disk reserve reached during preparation")
            while True:
                try:
                    result = results.get_nowait()
                except queue.Empty:
                    break
                if "error" in result:
                    if "MemoryError" in result["error"] or "OutOfMemory" in result["error"]:
                        raise MemoryError(result["error"])
                    raise RuntimeError(result["error"])
                if "shard" in result:
                    finished[result["worker"]] = result
                else:
                    total_rows += result["rows"]
                    notify(
                        "prepare",
                        f"{source}: {total_rows:,} rows prepared; {workers - len(finished)} workers",
                        rows=total_rows,
                    )
            if any(p.exitcode not in (None, 0) for p in processes):
                raise RuntimeError("Preparation worker exited unexpectedly")

        def submit(task):
            if task is not None:
                path, module, dataset, _mapped = task
                input_fingerprints[Path(path).name] = {
                    **file_fingerprint(path),
                    "module": module,
                    "dataset": dataset,
                }
            while True:
                drain()
                try:
                    tasks.put(task, timeout=0.2)
                    return
                except queue.Full:
                    pass

        try:
            for p in processes:
                p.start()
            batch_id = 0
            for ds in datasets:
                notify("prepare", f"Preparing {source}:{ds.dataset_name} with {workers} workers")
                raw = ds.raw_dataset
                iterator = (
                    raw.raw(
                        force_refresh=force_refresh,
                        source=source,
                        dataset=ds.dataset_name,
                        max_records=max_records,
                    )
                    if raw is not None
                    else ds.call()
                )
                handle = None
                dataset_rows = 0
                try:
                    for index, row in enumerate(islice(iterator, max_records)):
                        check_cancel(should_cancel)
                        dataset_rows += 1
                        if handle is None:
                            path = root / f"input-{batch_id}.pickle"
                            handle = path.open("wb")
                            count = 0
                        pickle.dump(
                            (index, row, canonical_json(row)),
                            handle,
                            protocol=pickle.HIGHEST_PROTOCOL,
                        )
                        count += 1
                        if count >= min(max_batch_records, 8192) or handle.tell() >= min(
                            raw_batch_bytes, 128 * 1024**2
                        ):
                            handle.close()
                            handle = None
                            submit(
                                (str(path), ds.qualified_module, ds.dataset_name, raw is not None)
                            )
                            batch_id += 1
                    if handle:
                        handle.close()
                        handle = None
                        submit((str(path), ds.qualified_module, ds.dataset_name, raw is not None))
                        batch_id += 1
                    if (
                        not dataset_rows
                        and raw is not None
                        and max_records != 0
                        and not getattr(raw, "allow_empty", False)
                    ):
                        # An empty dataset is almost always a failed download or
                        # a parser that no longer matches its input.
                        raise RuntimeError(
                            f"{source}:{ds.dataset_name} yielded no rows; check its "
                            "download, or declare the dataset allow_empty"
                        )
                finally:
                    if handle:
                        handle.close()
                    if hasattr(iterator, "close"):
                        iterator.close()
            for _ in processes:
                submit(None)
            while len(finished) < workers:
                drain()
                time.sleep(0.1)
            for p in processes:
                p.join()
            prepared = time.monotonic()
            check_cancel(should_cancel)
            notify(
                "finalize",
                f"Finalizing {source}: {workers} shards, {final_threads} SQL threads, "
                f"{final_memory} DuckDB memory",
            )
            writer = ParquetWriter(target, library_dir=library_dir, memory_limit=final_memory)
            writer.set_threads(final_threads)
            try:
                for i in sorted(finished):
                    drain()
                    writer.import_observation_shard(finished[i]["shard"])
                paths = [
                    str(Path(finished[i]["shard"]["directory"]) / "resolution_keys.parquet")
                    for i in sorted(finished)
                ]
                resolution = writer.resolution_summary([Path(path) for path in paths])
                resolution.update(
                    library_dir=str(library_dir) if library_dir else None,
                    libraries=finished[0]["resolution"].get("libraries", []),
                    lookup_metrics={},
                )
                for info in finished.values():
                    for key, value in info["resolution"].get("lookup_metrics", {}).items():
                        resolution["lookup_metrics"][key] = (
                            resolution["lookup_metrics"].get(key, 0) + value
                        )
                notify("finalize", f"Aggregating {source} observations and writing Parquets")
                outputs = writer.close()
            except BaseException:
                writer.abort()
                raise
            resolution["complex_composition"] = {
                "scope": "unique_precomposition_output_keys",
                "resolved": writer.metrics.get("composition_complexes", 0),
            }
            resolution_path = target / "resolution_stats.json"
            resolution_path.write_text(json.dumps(resolution, indent=2) + "\n")
            elapsed = time.monotonic() - started
            return dict(
                resource=source,
                input_fingerprints=input_fingerprints,
                files=outputs["files"],
                rows=outputs["rows"],
                resolution_stats_path=resolution_path,
                resolution_stats=resolution,
                elapsed_seconds=elapsed,
                phase_metrics={
                    "prepare_seconds": prepared - started,
                    "finalize_seconds": time.monotonic() - prepared,
                },
                batch_metrics={
                    key: (
                        sum(info["batch_metrics"][key] for info in finished.values())
                        if key == "flushes"
                        else max(info["batch_metrics"][key] for info in finished.values())
                    )
                    for key in ("flushes", "peak_records", "peak_relations", "peak_estimated_bytes")
                },
                writer_metrics=writer.metrics,
                duckdb_memory_limit=final_memory,
                duckdb_threads=final_threads,
                build_execution="parallel-shards-v1",
                batch_workers=workers,
                batch_limits={
                    "records": min(max_batch_records, 8192),
                    "raw_bytes": min(raw_batch_bytes, 128 * 1024**2),
                    "entities": entity_limit,
                    "relations": relation_limit,
                    "estimated_bytes": limits["bytes"],
                },
                payload_origins={
                    ds.dataset_name: "parsed_record"
                    if ds.raw_dataset is not None
                    else "mapped_record_fallback"
                    for ds in datasets
                },
            )
        finally:
            for p in processes:
                if p.is_alive():
                    p.terminate()
            for p in processes:
                if p.pid:
                    p.join(timeout=2)
                    if p.is_alive():
                        p.kill()
                        p.join()
            tasks.cancel_join_thread()
            tasks.close()
            results.close()
