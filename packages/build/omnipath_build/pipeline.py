"""High-level pipeline orchestrating inputs_v2 ingestion directly into serving Parquet files.

Discovers datasets, streams raw records through the Silver extractor, performs batch
entity resolution, appends observations, and merges once into the 3 required Parquet files.
"""

from __future__ import annotations

import concurrent.futures
import hashlib
import json
import logging
import tempfile
from datetime import datetime, timezone
from collections.abc import Sequence
from pathlib import Path
from typing import Any, Unpack

from .contracts import BuildConfig, BuildOptions, BuildResult, TwoPhaseConfig
from .cachedir_compat import patch_cachedir_opener
from .discovery import discover_datasets
from .progress import BuildCancelled, CancelCheck, ProgressCallback, check_cancel, emit
from omnipath_resolver import locate_library_dir
from omnipath_core.versioning import (
    BUILD_SCHEMA_VERSION,
    RESOURCE_FILES,
    SERVING_FILES,
    SERVING_SCHEMA_VERSION,
    validate_build_manifest,
    validate_source,
    validate_version,
)

logger = logging.getLogger(__name__)


def build_resource(
    source: str,
    *,
    version: str | None = None,
    output_dir: str | Path = "data",
    config: BuildConfig | None = None,
    **options: Unpack[BuildOptions],
) -> BuildResult:
    """Build an explicit, immutable resource version; publish only on success."""
    kwargs = BuildConfig(**{**(config.options() if config else {}), **options}).options()
    patch_cachedir_opener()
    import pyarrow.parquet as pq

    version = validate_version(version)
    source = validate_source(source.strip().lower())
    root = Path(output_dir).resolve()
    parent = root / "resources" / source
    parent.mkdir(parents=True, exist_ok=True)
    target = parent / version
    lock = parent / f".{version}.lock"
    from .locking import BuildLock

    with BuildLock(lock):
        if target.exists():
            raise FileExistsError(
                f"Resource {source}/{version} already exists; choose a new version"
            )
        with tempfile.TemporaryDirectory(prefix=".build-", dir=root) as staging:
            kwargs["library_dir"] = locate_library_dir(kwargs.get("library_dir"), root)
            result = _build_resource(source, version=version, output_dir=staging, **kwargs)
            if result["resource"] != source:
                raise ValueError("Discovered resource name differs from requested source")
            staged = Path(staging) / "resources" / source / version
            files = {}
            for name in RESOURCE_FILES + SERVING_FILES:
                path = staged / name
                rows = pq.read_metadata(path).num_rows
                digest = hashlib.sha256()
                with path.open("rb") as handle:
                    for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                        digest.update(chunk)
                files[name] = {
                    "rows": rows,
                    "sha256": digest.hexdigest(),
                    "size_bytes": path.stat().st_size,
                }
            from .provenance import runtime_provenance, reference_provenance

            manifest = {
                "schema_version": BUILD_SCHEMA_VERSION,
                "resource": source,
                "version": version,
                "provenance": runtime_provenance(),
                "reference": reference_provenance(kwargs["library_dir"]),
                "input_fingerprints": result.get("input_fingerprints", {}),
                "relation_taxon_policy": "asserted-consensus-else-participant-consensus-v1",
                "created_at": datetime.now(timezone.utc).isoformat(),
                "max_records": kwargs.get("max_records"),
                "datasets": kwargs.get("datasets"),
                "processing_policy": "molecular-form-observations-v3",
                "serving_schema_version": SERVING_SCHEMA_VERSION,
                "entity_resolution_policy": "gene-reference-product-identity-v3",
                "complex_resolution_policy": "composition-molecular-set-v2",
                "batch_limits": result.get("batch_limits"),
                "build_execution": result.get("build_execution"),
                "batch_workers": result.get("batch_workers", 0),
                "duckdb_memory_limit": result.get("duckdb_memory_limit"),
                "duckdb_threads": result.get("duckdb_threads"),
                "payload_origins": result.get("payload_origins"),
                "phase_metrics": result.get("phase_metrics"),
                "resource_metadata": result.get("resource_metadata"),
                "files": files,
            }
            validate_build_manifest(manifest)
            (staged / "build_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
            check_cancel(kwargs.get("should_cancel"))
            staged.rename(target)
            result["files"] = {
                name: target / Path(path).name for name, path in result["files"].items()
            }
            result["resolution_stats_path"] = target / Path(result["resolution_stats_path"]).name
            result["manifest_path"] = target / "build_manifest.json"
            return result


def _notify(
    on_progress: ProgressCallback | None,
    *,
    prints: bool,
    **event: Any,
) -> None:
    emit(on_progress, **event)
    message = event.get("message")
    if prints and message:
        print(message, flush=True)


def _build_resource(
    source: str,
    *,
    version: str | None = None,
    output_dir: str | Path = "data",
    max_records: int | None = None,
    batch_size: int = 50000,
    batch_workers: int | None = None,
    resource_ram_bytes: int | None = None,
    min_free_disk_bytes: int = 20 * 1024**3,
    raw_batch_bytes: int | None = None,
    duckdb_memory_limit: str | None = None,
    max_batch_records: int = 50000,
    max_batch_relations: int = 50000,
    max_batch_bytes: int = 256 * 1024 * 1024,
    force_refresh: bool = False,
    progress: bool = True,
    inputs_package: str = "pypath.inputs_v2",
    datasets: Sequence[str] | None = None,
    cache_dir: str | Path | None = None,
    on_progress: ProgressCallback | None = None,
    should_cancel: CancelCheck | None = None,
    library_dir: str | Path | None = None,
) -> dict[str, Any]:
    """Ingest, canonicalize, and write Parquet tables for a single resource using streaming batches."""
    from .duckdb_config import build_threads

    if batch_workers is None:
        import os

        available = (
            len(os.sched_getaffinity(0))
            if hasattr(os, "sched_getaffinity")
            else (os.cpu_count() or 1)
        )
        batch_workers = max(1, min(build_threads(), available - 1))
    if batch_workers < 1:
        raise ValueError("At least one preparation worker is required; use batch_workers=1")
    if raw_batch_bytes is not None and raw_batch_bytes < 1:
        raise ValueError("Raw batch byte limit must be positive")
    if min(batch_size, max_batch_records, max_batch_relations, max_batch_bytes) < 1:
        raise ValueError("All batch limits must be positive")

    source_stage = f"build:{source.lower()}"

    _notify(
        on_progress,
        prints=progress,
        pipeline="resource",
        stage="discover",
        status="running",
        message=f"Discovering datasets for {source}",
    )
    check_cancel(should_cancel)
    source_slug, discovered_list, _config = discover_datasets(
        source=source,
        inputs_package=inputs_package,
        datasets=datasets,
        cache_dir=cache_dir,
    )
    _notify(
        on_progress,
        prints=progress,
        pipeline="resource",
        stage="discover",
        status="done",
        current=len(discovered_list),
        message=f"Found {len(discovered_list)} dataset(s) for {source_slug}",
    )

    version = validate_version(version)

    _notify(
        on_progress,
        prints=progress,
        pipeline="resource",
        stage=source_stage,
        status="running",
        current=0,
        total=len(discovered_list),
        message=f"Building {source_slug} ({len(discovered_list)} dataset(s), batch_size={batch_size:,})",
        counts={"entities": 0, "relations": 0, "payloads": 0, "rows": 0},
    )

    target_dir = Path(output_dir).resolve() / "resources" / source_slug / version
    from .two_phase import build_two_phase
    from .duckdb_config import build_memory_limit, finalize_budget

    ram = resource_ram_bytes or max(1024**3, batch_workers * 1024**3)
    if resource_ram_bytes:
        final_memory, final_threads = finalize_budget(resource_ram_bytes, build_threads())
    else:
        final_memory, final_threads = build_memory_limit(duckdb_memory_limit), None
    result = build_two_phase(
        source_slug,
        discovered_list,
        target_dir,
        config=TwoPhaseConfig(
            workers=batch_workers,
            cpus=build_threads(),
            ram_bytes=ram,
            final_memory=final_memory,
            library_dir=library_dir,
            max_records=max_records,
            force_refresh=force_refresh,
            max_batch_records=max_batch_records,
            raw_batch_bytes=raw_batch_bytes or max_batch_bytes,
            min_free_disk_bytes=min_free_disk_bytes,
            entity_limit=min(batch_size, 5000),
            relation_limit=min(max_batch_relations, 5000),
            byte_limit=max_batch_bytes,
            final_threads=final_threads,
        ),
        on_progress=on_progress,
        should_cancel=should_cancel,
    )
    from .source_metadata import snapshot_source_metadata
    from .discovery import setup_pypath_cache

    result["resource_metadata"] = snapshot_source_metadata(
        source_slug, discovered_list, setup_pypath_cache(cache_dir)
    )
    result["version"] = version
    return result


class BuildBatchError(RuntimeError):
    """Completed batch with unsuccessful resources; successful files remain valid."""

    def __init__(self, results: dict[str, dict[str, Any]]):
        self.results = results
        failed = [key for key, value in results.items() if value.get("status") == "failed"]
        super().__init__(
            "Resource builds failed (remaining sources were attempted): " + ", ".join(failed)
        )


def build_all(
    sources: Sequence[str],
    *,
    version: str | None = None,
    versions: dict[str, str] | None = None,
    output_dir: str | Path = "data",
    max_records: int | None = None,
    batch_size: int = 50000,
    duckdb_memory_limit: str | None = None,
    force_refresh: bool = False,
    progress: bool = True,
    parallel: int = 1,
    on_progress: ProgressCallback | None = None,
    should_cancel: CancelCheck | None = None,
) -> dict[str, dict[str, Any]]:
    """Build all specified sources in sequence or in parallel."""
    results: dict[str, dict[str, Any]] = {}
    if not sources:
        return results
    sources = list(dict.fromkeys(validate_source(src.strip().lower()) for src in sources))
    selected_versions = {
        src: validate_version((versions or {}).get(src, version)) for src in sources
    }
    for src, ver in selected_versions.items():
        if (Path(output_dir) / "resources" / src / ver).exists():
            raise FileExistsError(f"Resource {src}/{ver} already exists; choose a new version")

    def _build_one(src: str) -> dict[str, Any]:
        check_cancel(should_cancel)
        try:
            return build_resource(
                source=src,
                version=selected_versions[src],
                output_dir=output_dir,
                max_records=max_records,
                batch_size=batch_size,
                duckdb_memory_limit=duckdb_memory_limit,
                force_refresh=force_refresh,
                progress=progress,
                on_progress=on_progress,
                should_cancel=should_cancel,
            )
        except BuildCancelled:
            raise
        except Exception as exc:
            logger.exception("Skipping failed resource %s; continuing batch", src)
            emit(
                on_progress,
                stage=f"build:{src}",
                status="failed",
                message=f"{src}: {type(exc).__name__}: {exc}",
            )
            return {
                "resource": src,
                "version": selected_versions[src],
                "status": "failed",
                "error": f"{type(exc).__name__}: {exc}",
            }

    if parallel <= 1 or len(sources) <= 1:
        for src in sources:
            res = _build_one(src)
            results[f"{res['resource']}/{res['version']}"] = res
    else:
        workers = min(int(parallel), len(sources))
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
            futures = [executor.submit(_build_one, src) for src in sources]
            for future in concurrent.futures.as_completed(futures):
                check_cancel(should_cancel)
                res = future.result()
                results[f"{res['resource']}/{res['version']}"] = res

    if any(result.get("status") == "failed" for result in results.values()):
        raise BuildBatchError(results)
    return results
