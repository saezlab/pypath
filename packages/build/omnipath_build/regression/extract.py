"""Extract pre-resolution observations from raw source extraction.

For each resource, every dataset is run through the normal pypath extraction
(the same raw iterator and mapper the build uses), records go through
``SilverExtractor`` and each extracted entity is turned into ``queries`` and
``votes`` by ``votes_for`` and ``observation_bundle``, exactly as
``LibraryMatcher.targets`` does at build time. Nothing published by a build is
read.
"""

from __future__ import annotations

import concurrent.futures
import json
import multiprocessing
import os
import time
import traceback
from collections import Counter
from collections.abc import Iterable, Iterator, Mapping
from itertools import islice
from pathlib import Path
from typing import Any

from .store import (
    ObservationWriter,
    observation_fingerprint,
    resource_dir,
)

Observation = tuple[str, dict[str, Any], list[dict[str, Any]]]


def log(event: str, **kw: Any) -> None:
    print(
        json.dumps({"time": time.strftime("%FT%TZ", time.gmtime()), "event": event, **kw}),
        flush=True,
    )


def entity_observations(entities: Mapping[str, Any]) -> tuple[list[Observation], int]:
    """(fingerprint, query, vote rows) per extracted entity with a reference library.

    Entities without a library (complexes, CV terms, ...) never reach the
    resolver and are only counted.
    """
    from omnipath_resolver.canonical.match import votes_for
    from omnipath_resolver.canonical.policy import get_policy
    from omnipath_resolver.observations import observation_bundle

    out: list[Observation] = []
    skipped = 0
    for key, obs in entities.items():
        policy = get_policy(obs.entity_type)
        if policy.library is None:
            skipped += 1
            continue
        normalized, observed = votes_for(obs, policy)
        query, rows = observation_bundle(key, obs, normalized, observed, policy.library)
        query["taxon"] = str(getattr(obs, "taxon", None) or "")
        for row in rows:
            row.setdefault("gene_only", False)
        out.append((observation_fingerprint(query["target"], rows), query, rows))
    return out, skipped


def process_records(
    source: str,
    dataset: str,
    records: Iterable[tuple[int, Any]],
    mapper: Any = None,
) -> tuple[list[Observation], int]:
    """Run raw records (index, row) through the mapper and the silver extractor."""
    from omnipath_resolver.observations import flush_lipid_cache

    from ..silver import SilverExtractor

    extractor = SilverExtractor(source, dataset)
    for index, raw in records:
        record = mapper(raw) if mapper else raw
        if record is not None:
            extractor.process_record(record, raw, f"{dataset}:{index}", index, payload_json="")
    result = entity_observations(extractor.entities)
    flush_lipid_cache()
    return result


def process_batch(task: tuple[str, str, str, bool, int, list[Any]]) -> dict[str, Any]:
    """Worker entry point: one batch of raw rows of one dataset."""
    source, module, dataset, mapped, start, rows = task
    from ..cachedir_compat import patch_cachedir_opener
    from ..two_phase import load_mapper

    patch_cachedir_opener()
    mapper = load_mapper(module, dataset) if mapped else None
    try:
        observations, skipped = process_records(source, dataset, enumerate(rows, start), mapper)
    except Exception as exc:
        raise RuntimeError(
            f"{source}:{dataset} rows {start}-{start + len(rows) - 1} ({module}) failed: "
            f"{type(exc).__name__}: {exc}"
        ) from exc
    return dict(rows=len(rows), observations=observations, skipped=skipped)


def _row_batches(
    iterator: Iterator[Any], max_records: int | None, size: int
) -> Iterator[tuple[int, list[Any]]]:
    start = 0
    limited = islice(iterator, max_records)
    while True:
        rows = list(islice(limited, size))
        if not rows:
            return
        yield start, rows
        start += len(rows)


class _Collector:
    """Global dedupe across datasets, with per-dataset statistics."""

    def __init__(self, writer: ObservationWriter):
        self.writer = writer
        self.seen: set[str] = set()
        self.dataset: dict[str, dict[str, Any]] = {}

    def start(self, dataset: str) -> dict[str, Any]:
        stats = self.dataset[dataset] = dict(
            status="running",
            records=0,
            observations=0,
            new_observations=0,
            skipped_no_library=0,
            seconds=0.0,
            error=None,
        )
        self._occurrences: Counter[str] = Counter()
        self._current = dataset
        return stats

    def absorb(self, result: dict[str, Any]) -> None:
        stats = self.dataset[self._current]
        stats["records"] += result["rows"]
        stats["skipped_no_library"] += result["skipped"]
        for fingerprint, query, rows in result["observations"]:
            stats["observations"] += 1
            self._occurrences[fingerprint] += 1
            if fingerprint not in self.seen:
                self.seen.add(fingerprint)
                stats["new_observations"] += 1
                self.writer.add(fingerprint, query, rows)
        if len(self._occurrences) >= 500_000:
            self.flush_occurrences()

    def flush_occurrences(self) -> None:
        self.writer.add_occurrences(self._current, dict(self._occurrences))
        self._occurrences = Counter()


def extract_resource(
    source: str,
    output: str | Path,
    *,
    max_records: int | None = None,
    workers: int = 2,
    batch_records: int = 2000,
    cache_dir: str | Path | None = None,
    datasets: list[str] | None = None,
    memory_limit: str = "2GB",
) -> dict[str, Any]:
    """Extract one resource into ``<output>/<source>``; never raises for source errors."""
    from ..cachedir_compat import patch_cachedir_opener
    from ..discovery import discover_datasets

    patch_cachedir_opener()
    started = time.monotonic()
    directory = resource_dir(output, source)
    info: dict[str, Any] = dict(
        resource=source,
        status="failed",
        max_records=max_records,
        workers=workers,
        batch_records=batch_records,
        pypath_cache=str(cache_dir or os.environ.get("PYPATH_DOWNLOAD_DATADIR", "")),
        started=time.strftime("%FT%TZ", time.gmtime()),
        datasets={},
        error=None,
    )
    # A stale partial extraction must never be mistaken for a finished one.
    for stale in ("queries", "votes", "occurrences"):
        (directory / f"{stale}.parquet").unlink(missing_ok=True)
    directory.mkdir(parents=True, exist_ok=True)
    # A process killed mid-run (e.g. out of memory) leaves this marker behind.
    info["status"] = "running"
    (directory / "extract.json").write_text(json.dumps(info, indent=2, default=str) + "\n")
    info["status"] = "failed"
    writer = collector = None
    try:
        _, discovered, _ = discover_datasets(source, datasets=datasets, cache_dir=cache_dir)
        writer = ObservationWriter(directory)
        collector = _Collector(writer)
        pool = (
            concurrent.futures.ProcessPoolExecutor(
                max_workers=workers, mp_context=multiprocessing.get_context("spawn")
            )
            if workers > 1
            else None
        )
        try:
            for ds in discovered:
                name = ds.dataset_name
                stats = collector.start(name)
                t0 = time.monotonic()
                log("dataset_start", resource=source, dataset=name)
                try:
                    _run_dataset(source, ds, collector, pool, workers, max_records, batch_records)
                    stats["status"] = "ok"
                except Exception as exc:
                    stats["status"] = "failed"
                    stats["error"] = f"{type(exc).__name__}: {exc}"
                    stats["traceback"] = traceback.format_exc()[-4000:]
                    log("dataset_failed", resource=source, dataset=name, error=stats["error"])
                finally:
                    collector.flush_occurrences()
                    stats["seconds"] = round(time.monotonic() - t0, 2)
                log(
                    "dataset_done",
                    resource=source,
                    dataset=name,
                    **{k: v for k, v in stats.items() if k != "traceback"},
                )
        finally:
            if pool is not None:
                pool.shutdown(wait=True, cancel_futures=True)
        counts = writer.finalize(memory_limit=memory_limit)
        info["rows"] = counts
        info["datasets"] = collector.dataset
        failed = [n for n, s in collector.dataset.items() if s["status"] != "ok"]
        info["status"] = "ok" if not failed else ("partial" if counts["queries"] else "failed")
        if failed:
            info["error"] = "datasets failed: " + ", ".join(failed)
        info["observations_by_library"] = _library_counts(directory / "queries.parquet")
    except Exception as exc:
        info["error"] = f"{type(exc).__name__}: {exc}"
        info["traceback"] = traceback.format_exc()[-4000:]
        if collector is not None:
            info["datasets"] = collector.dataset
    info["seconds"] = round(time.monotonic() - started, 2)
    (directory / "extract.json").write_text(json.dumps(info, indent=2, default=str) + "\n")
    log(
        "resource_done",
        resource=source,
        status=info["status"],
        seconds=info["seconds"],
        by_library=info.get("observations_by_library"),
        error=info["error"],
    )
    return info


def _library_counts(queries_path: Path) -> dict[str, int]:
    import duckdb

    con = duckdb.connect()
    rows = con.execute(
        f"SELECT library, count(*) FROM read_parquet('{queries_path}') GROUP BY library"
    ).fetchall()
    con.close()
    return {library: n for library, n in sorted(rows)}


def _run_dataset(
    source: str,
    ds: Any,
    collector: _Collector,
    pool: concurrent.futures.Executor | None,
    workers: int,
    max_records: int | None,
    batch_records: int,
) -> None:
    from ..two_phase import load_mapper

    raw = ds.raw_dataset
    iterator = (
        raw.raw(
            force_refresh=False,
            source=source,
            dataset=ds.dataset_name,
            max_records=max_records,
        )
        if raw is not None
        else ds.call()
    )
    mapped = raw is not None
    pending: set[concurrent.futures.Future] = set()
    failure: BaseException | None = None

    def drain(limit: int) -> None:
        nonlocal pending, failure
        while len(pending) > limit:
            done, pending = concurrent.futures.wait(
                pending, return_when=concurrent.futures.FIRST_COMPLETED
            )
            for future in done:
                try:
                    collector.absorb(future.result())
                except BaseException as exc:  # first failure wins; the rest drain
                    failure = failure or exc

    try:
        mapper = (
            load_mapper(ds.qualified_module, ds.dataset_name) if (mapped and pool is None) else None
        )
        for start, rows in _row_batches(iterator, max_records, batch_records):
            if failure is not None:
                break
            if pool is None:
                collector.absorb(
                    _inline(source, ds.dataset_name, start, rows, mapper, ds.qualified_module)
                )
            else:
                pending.add(
                    pool.submit(
                        process_batch,
                        (source, ds.qualified_module, ds.dataset_name, mapped, start, rows),
                    )
                )
                drain(workers * 2)
        drain(0)
    finally:
        if hasattr(iterator, "close"):
            iterator.close()
        for future in pending:
            future.cancel()
        drain(0)
    if failure is not None:
        raise failure


def _inline(source, dataset, start, rows, mapper, module):
    try:
        observations, skipped = process_records(source, dataset, enumerate(rows, start), mapper)
    except Exception as exc:
        raise RuntimeError(
            f"{source}:{dataset} rows {start}-{start + len(rows) - 1} ({module}) failed: "
            f"{type(exc).__name__}: {exc}"
        ) from exc
    return dict(rows=len(rows), observations=observations, skipped=skipped)


def extract_resources(
    resources: list[str],
    output: str | Path,
    *,
    skip_existing: bool = False,
    **options: Any,
) -> dict[str, dict[str, Any]]:
    """Extract several resources one after another; a failure never stops the rest."""
    results = {}
    os.environ.setdefault("OMNIPATH_GOSLIN_CACHE", str(Path(output).resolve() / ".goslin-cache"))
    for source in resources:
        marker = resource_dir(output, source) / "extract.json"
        if skip_existing and marker.is_file():
            previous = json.loads(marker.read_text())
            if previous.get("status") == "ok" and (
                previous.get("max_records") == options.get("max_records")
            ):
                results[source] = previous
                log("resource_skipped", resource=source)
                continue
        results[source] = extract_resource(source, output, **options)
    return results
