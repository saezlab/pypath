"""Resolve stored observations with a runtime and persist results plus timings."""

from __future__ import annotations

import json
import time
from collections import Counter
from pathlib import Path
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq

from .store import (
    RESULT_SCHEMA,
    iter_resolve_batches,
    list_resources,
    resource_dir,
)

IDENTITY_FORMAT = "omnipath-identity-v2"


def log(event: str, **kw: Any) -> None:
    print(
        json.dumps({"time": time.strftime("%FT%TZ", time.gmtime()), "event": event, **kw}),
        flush=True,
    )


def raise_open_file_limit() -> None:
    """Every hub index opens its own LMDB shards; lift the soft fd limit."""
    try:
        import resource

        soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
        wanted = hard if hard != resource.RLIM_INFINITY else 1 << 20
        if soft != resource.RLIM_INFINITY and soft < wanted:
            resource.setrlimit(resource.RLIMIT_NOFILE, (wanted, hard))
    except (ImportError, ValueError, OSError):
        pass


def load_runtime(path: str | Path, *, cache_dir: str | Path | None = None) -> tuple[Any, dict]:
    """Open the identity runtime of an ``omnipath-identity-v2`` directory."""
    from omnipath_resolver.identity_runtime import IdentityRuntime

    path = Path(path)
    raise_open_file_limit()
    manifest = json.loads((path / "manifest.json").read_text())
    if manifest.get("format") != IDENTITY_FORMAT:
        raise ValueError(f"{path} is not an identity directory ({IDENTITY_FORMAT})")
    return IdentityRuntime(path, cache_dir=cache_dir), manifest


def result_rows(resolved: dict[str, Any]) -> list[dict[str, Any]]:
    """Flatten ``runtime.resolve`` output into RESULT_SCHEMA rows."""
    records = resolved.get("records", {})
    rows = []
    for r in resolved["results"]:
        entities = list(r.get("entities") or [])
        rows.append(
            dict(
                input_id=r["input_id"],
                outcome=str(r["outcome"]),
                entities=entities,
                gene_mapping_status=r.get("gene_mapping_status"),
                gene_candidates=[str(g) for g in (r.get("gene_candidates") or [])],
                protein_entity_id=r.get("protein_entity_id") or None,
                candidate_count=int(r.get("candidate_count") or 0),
                entity_labels=[str((records.get(e) or {}).get("label") or "") for e in entities],
            )
        )
    return rows


def resolve_resource(
    runtime: Any,
    observations: str | Path,
    results: str | Path,
    resource: str,
    *,
    batch_size: int = 5000,
    sample: int | None = None,
) -> dict[str, Any]:
    """Resolve one resource's stored observations in batches; write results Parquet."""
    source = resource_dir(observations, resource)
    target = resource_dir(results, resource)
    target.mkdir(parents=True, exist_ok=True)
    tmp = target / ".results.parquet.tmp"
    batches: list[dict[str, Any]] = []
    outcomes: Counter[tuple[str, str]] = Counter()
    wall_start = time.perf_counter()
    totals = Counter()
    with pq.ParquetWriter(tmp, RESULT_SCHEMA, compression="zstd") as writer:
        for number, (queries, votes) in enumerate(
            iter_resolve_batches(source, batch_size, sample=sample)
        ):
            started = time.perf_counter()
            resolved, metrics = runtime.resolve(queries, votes)
            wall = time.perf_counter() - started
            rows = result_rows(resolved)
            if len(rows) != len(queries):
                raise RuntimeError(
                    f"{resource}: runtime returned {len(rows)} results for {len(queries)} queries"
                )
            writer.write_table(pa.Table.from_pylist(rows, schema=RESULT_SCHEMA))
            library = {q["input_id"]: q["library"] for q in queries}
            for row in rows:
                outcomes[(library[row["input_id"]], row["outcome"])] += 1
            batch = dict(
                batch=number,
                queries=len(queries),
                votes=len(votes),
                lookup_seconds=metrics.get("lookup_seconds", 0.0),
                decision_seconds=metrics.get("decision_seconds", 0.0),
                record_seconds=metrics.get("entity_fetch_seconds", 0.0),
                total_seconds=metrics.get("total_seconds", wall),
                wall_seconds=wall,
                unique_keys=metrics.get("unique_keys"),
                candidate_records=metrics.get("candidate_records"),
                entity_records=metrics.get("entity_records"),
            )
            batches.append(batch)
            totals.update(
                queries=len(queries),
                votes=len(votes),
                lookup_seconds=batch["lookup_seconds"],
                decision_seconds=batch["decision_seconds"],
                record_seconds=batch["record_seconds"],
                wall_seconds=wall,
            )
    tmp.replace(target / "results.parquet")
    elapsed = time.perf_counter() - wall_start
    summary = dict(
        resource=resource,
        observations=totals["queries"],
        votes=totals["votes"],
        lookup_seconds=totals["lookup_seconds"],
        decision_seconds=totals["decision_seconds"],
        record_seconds=totals["record_seconds"],
        resolve_seconds=totals["wall_seconds"],
        elapsed_seconds=elapsed,
        observations_per_second=(totals["queries"] / totals["wall_seconds"])
        if totals["wall_seconds"]
        else None,
        outcomes={
            library: {o: n for (lib, o), n in sorted(outcomes.items()) if lib == library}
            for library in sorted({lib for lib, _ in outcomes})
        },
        batches=batches,
    )
    log(
        "resource_resolved",
        resource=resource,
        observations=summary["observations"],
        seconds=round(elapsed, 2),
        observations_per_second=summary["observations_per_second"],
        outcomes=summary["outcomes"],
    )
    return summary


def resolve_all(
    runtime_path: str | Path,
    observations: str | Path,
    results: str | Path,
    *,
    resources: list[str] | None = None,
    batch_size: int = 5000,
    cache_dir: str | Path | None = None,
    force: bool = False,
    runtime: Any = None,
    sample: int | None = None,
) -> dict[str, Any]:
    """Resolve every selected resource; metrics.json accumulates across runs."""
    results = Path(results)
    results.mkdir(parents=True, exist_ok=True)
    available = list_resources(observations)
    selected = resources or available
    missing = sorted(set(selected) - set(available))
    if missing:
        raise FileNotFoundError(f"No extracted observations for: {', '.join(missing)}")
    metrics_path = results / "metrics.json"
    metrics = json.loads(metrics_path.read_text()) if metrics_path.is_file() else {}
    started = time.perf_counter()
    owned = runtime is None
    manifest: dict[str, Any] = {}
    if owned:
        runtime, manifest = load_runtime(runtime_path, cache_dir=cache_dir)
    load_seconds = time.perf_counter() - started
    metrics.update(
        runtime=dict(
            path=str(runtime_path),
            class_name=type(runtime).__name__,
            format=manifest.get("format"),
            fingerprint=manifest.get("fingerprint") or manifest.get("reference_fingerprint"),
            load_seconds=load_seconds,
        ),
        observations=str(observations),
        batch_size=batch_size,
        sample_per_resource_and_library=sample or None,
    )
    per_resource = metrics.setdefault("resources", {})
    try:
        for resource in selected:
            if not force and (resource_dir(results, resource) / "results.parquet").is_file():
                log("resource_skipped", resource=resource)
                continue
            try:
                per_resource[resource] = resolve_resource(
                    runtime, observations, results, resource, batch_size=batch_size, sample=sample
                )
            except Exception as exc:
                import traceback

                per_resource[resource] = dict(
                    resource=resource,
                    error=f"{type(exc).__name__}: {exc}",
                    traceback=traceback.format_exc()[-4000:],
                )
                log("resource_failed", resource=resource, error=per_resource[resource]["error"])
            finally:
                _write_metrics(metrics_path, metrics)
    finally:
        if owned and hasattr(runtime, "close"):
            runtime.close()
    return metrics


def _write_metrics(path: Path, metrics: dict[str, Any]) -> None:
    done = [r for r in metrics.get("resources", {}).values() if "observations" in r]
    wall = sum(r["resolve_seconds"] for r in done)
    observations = sum(r["observations"] for r in done)
    metrics["totals"] = dict(
        resources=len(done),
        observations=observations,
        votes=sum(r["votes"] for r in done),
        lookup_seconds=sum(r["lookup_seconds"] for r in done),
        decision_seconds=sum(r["decision_seconds"] for r in done),
        record_seconds=sum(r["record_seconds"] for r in done),
        resolve_seconds=wall,
        observations_per_second=(observations / wall) if wall else None,
    )
    tmp = path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(metrics, indent=2, default=str) + "\n")
    tmp.replace(path)
