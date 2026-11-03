"""Public build configuration and portable observation shard contracts."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, NotRequired, TypedDict
from collections.abc import Sequence
from .progress import CancelCheck, ProgressCallback


@dataclass(frozen=True)
class BuildConfig:
    max_records: int | None = None
    batch_size: int = 50000
    batch_workers: int | None = None
    resource_ram_bytes: int | None = None
    min_free_disk_bytes: int = 20 * 1024**3
    raw_batch_bytes: int | None = None
    duckdb_memory_limit: str | None = None
    max_batch_records: int = 50000
    max_batch_relations: int = 50000
    max_batch_bytes: int = 256 * 1024 * 1024
    force_refresh: bool = False
    progress: bool = True
    inputs_package: str = "pypath.inputs_v2"
    datasets: Sequence[str] | None = None
    cache_dir: str | Path | None = None
    on_progress: ProgressCallback | None = None
    should_cancel: CancelCheck | None = None
    library_dir: str | Path | None = None

    def options(self) -> dict[str, Any]:
        return {name: getattr(self, name) for name in self.__dataclass_fields__}


class ObservationShard(TypedDict):
    directory: str
    tables: list[str]
    events: int
    payloads: int
    metrics: dict[str, int]


class WorkerResult(TypedDict):
    worker: int
    shard: ObservationShard
    resolution: dict[str, Any]
    seconds: float
    batch_metrics: dict[str, int]


class BuildResult(TypedDict):
    resource: str
    files: dict[str, Path]  # table name -> Parquet path (resource and serving tables)
    rows: dict[str, int]
    resolution_stats_path: Path
    resolution_stats: dict[str, Any]
    elapsed_seconds: float
    phase_metrics: dict[str, float]
    batch_metrics: dict[str, int]
    writer_metrics: dict[str, int]
    duckdb_memory_limit: str
    duckdb_threads: int
    build_execution: str
    batch_workers: int
    batch_limits: dict[str, int]
    payload_origins: dict[str, str]
    input_fingerprints: dict[str, dict[str, Any]]
    version: NotRequired[str]
    resource_metadata: NotRequired[dict[str, Any]]
    manifest_path: NotRequired[Path]


class BuildOptions(TypedDict, total=False):
    """Typed keyword adapter for callers migrating to BuildConfig."""

    max_records: int | None
    batch_size: int
    batch_workers: int | None
    resource_ram_bytes: int | None
    min_free_disk_bytes: int
    raw_batch_bytes: int | None
    duckdb_memory_limit: str | None
    max_batch_records: int
    max_batch_relations: int
    max_batch_bytes: int
    force_refresh: bool
    progress: bool
    inputs_package: str
    datasets: Sequence[str] | None
    cache_dir: str | Path | None
    on_progress: ProgressCallback | None
    should_cancel: CancelCheck | None
    library_dir: str | Path | None


@dataclass(frozen=True)
class TwoPhaseConfig:
    workers: int
    cpus: int
    ram_bytes: int
    final_memory: str
    library_dir: str | Path | None
    max_records: int | None
    force_refresh: bool
    max_batch_records: int
    raw_batch_bytes: int
    min_free_disk_bytes: int = 20 * 1024**3
    entity_limit: int = 5000
    relation_limit: int = 5000
    byte_limit: int = 64 * 1024**2
    # Finalizer threads; None uses every CPU (development runs without a hard memory limit).
    final_threads: int | None = None
