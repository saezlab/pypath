"""Admin job runner: inspect resolver/resources and trigger long-running builds."""

from __future__ import annotations

import json
import logging
import os
import re
import shutil
import threading
import time
import traceback
from collections.abc import Callable, Iterator
from pathlib import Path
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq

from omnipath_api.jobs.build_ops import BuildOps, BuildUnavailable, _import_build
from omnipath_api.jobs.store import JobStore
from omnipath_api.jobs.logs import read_log_tail
from omnipath_api.queries.store.inventory import ReleaseStore

from omnipath_api.jobs.runner import (
    Job,
    make_stage as _make_stage,
    set_stage_status as _set_stage_status,
)

logger = logging.getLogger(__name__)

TERMINAL_STATUSES = frozenset({"done", "failed", "cancelled"})
_SAFE_LOG_NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
_MAX_LOG_BYTES = 2_000_000
FALLBACK_HUB_NAMES: tuple[str, ...] = (
    "entrez",
    "uniprot",
    "ensg",
    "chebi",
    "pubchem",
    "refmet",
    "chembl",
    "hmdb",
    "lipidmaps",
    "swisslipids",
    "bigg",
    "metanetx",
    "mirbase",
    "taxon_species",
)
DICTIONARY_NAMES: tuple[str, ...] = ("id_type", "backend", "organism")
LIBRARY_NAMES: tuple[str, ...] = ("gene_protein", "chemical")


def hub_names() -> tuple[str, ...]:
    build = _import_build()
    if build is None:
        return FALLBACK_HUB_NAMES
    try:
        from omnipath_build.hubs.export import HUB_NAMES

        return tuple(HUB_NAMES)
    except ModuleNotFoundError as exc:
        if exc.name and (exc.name == "pypath" or exc.name.startswith("pypath.")):
            logger.info("Build parser dependency unavailable; using known hub names")
            return FALLBACK_HUB_NAMES
        raise


def _parquet_info(path: Path, name: str | None = None) -> dict[str, Any]:
    info: dict[str, Any] = {
        "name": name or path.stem,
        "exists": path.is_file(),
        "path": str(path),
        "rows": 0,
        "size_bytes": 0,
        "mtime": None,
    }
    if not path.is_file():
        return info
    stat = path.stat()
    info["size_bytes"] = stat.st_size
    info["mtime"] = stat.st_mtime
    try:
        info["rows"] = int(pq.read_metadata(path).num_rows)
    except (OSError, pa.ArrowException) as exc:
        logger.exception("Cannot inspect Parquet artifact %s", path)
        info["error"] = str(exc)
    return info


def _log_info(path: Path) -> dict[str, Any]:
    info: dict[str, Any] = {
        "exists": path.is_file(),
        "path": str(path),
        "size_bytes": 0,
        "mtime": None,
    }
    if not path.is_file():
        return info
    stat = path.stat()
    info["size_bytes"] = stat.st_size
    info["mtime"] = stat.st_mtime
    return info


def resolve_hubs_dir(data_root: Path) -> Path:
    env = os.environ.get("OMNIPATH_RESOLVER_DATA_DIR")
    if env:
        root = Path(env)
        hubs = root / "hubs"
        if hubs.is_dir():
            return hubs.resolve()
        return root.resolve()
    return (data_root / "reference" / "hubs").resolve()


class AdminService:
    def __init__(self, data_root: str | Path, ops: BuildOps | None = None) -> None:
        self.data_root = Path(data_root).resolve()
        self.ops = ops or BuildOps()
        self._injected_ops = ops is not None
        self.job_store = JobStore(self.data_root)
        self._lock = threading.Lock()
        self._jobs: dict[str, Job] = {}
        self._current_id: str | None = None
        self._sources_cache: tuple[float, list[dict[str, Any]]] | None = None
        self._cv = threading.Condition(self._lock)
        self._load_persisted_jobs()

    def hubs_dir(self) -> Path:
        return resolve_hubs_dir(self.data_root)

    def library_dir(self) -> Path:
        env = os.environ.get("OMNIPATH_LIBRARY_DIR")
        if env:
            return Path(env).resolve()
        return self.hubs_dir().parent / "library"

    def build_capabilities(self, *, use_worker=True):
        if use_worker:
            persisted = self.job_store.capabilities()
            if persisted is not None:
                return persisted
        return {
            "build_available": self.ops.available(),
            "hubs_available": getattr(self.ops, "hubs_available", self.ops.available)(),
            "pipeline_available": getattr(self.ops, "pipeline_available", self.ops.available)(),
        }

    def inspect(self) -> dict[str, Any]:
        hubs_dir = self.hubs_dir()
        names = hub_names()
        hub_files = [
            {
                **_parquet_info(hubs_dir / f"{name}.parquet", name),
                "log": _log_info(self._hub_log_path(name)),
            }
            for name in names
        ]
        dictionaries = [
            _parquet_info(hubs_dir / f"{name}.parquet", name) for name in DICTIONARY_NAMES
        ]
        manifest_path = hubs_dir / "manifest.json"
        manifest: dict[str, Any] | None = None
        if manifest_path.is_file():
            try:
                manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            except (OSError, ValueError):
                logger.exception("Cannot read hub manifest %s", manifest_path)
                manifest = None
        library_dir = self.library_dir().resolve()
        library_manifest: dict[str, Any] | None = None
        if (library_dir / "manifest.json").is_file():
            try:
                library_manifest = json.loads(
                    (library_dir / "manifest.json").read_text(encoding="utf-8")
                )
            except (OSError, ValueError):
                logger.exception("Cannot read library manifest %s", library_dir)
                library_manifest = None
        library_ready = bool(
            library_manifest and library_manifest.get("format") == "omnipath-identity-v2"
        )
        built = self._built_resources()
        self._refresh_jobs()
        with self._lock:
            current = None
            if self._current_id and self._current_id in self._jobs:
                current = self._jobs[self._current_id].snapshot()
            jobs = [item.snapshot() for item in self._recent_jobs()]
        return {
            "data_root": str(self.data_root),
            **self.build_capabilities(),
            "resolver": {
                "ready": library_ready,
                "hubs_dir": str(hubs_dir),
                "library_dir": str(library_dir),
                "library_manifest": library_manifest,
                "hubs": hub_files,
                "dictionaries": dictionaries,
                "manifest": manifest,
            },
            "resources": {
                "built": built,
                "logs": self._resource_log_map(),
            },
            "job": current,
            "jobs": jobs,
        }

    def available_sources(self) -> list[dict[str, Any]]:
        return self._available_sources()

    def current_job(self) -> dict[str, Any] | None:
        self._refresh_jobs()
        with self._lock:
            if self._current_id and self._current_id in self._jobs:
                return self._jobs[self._current_id].snapshot()
            return None

    def get_job(self, job_id: str) -> dict[str, Any] | None:
        self._refresh_jobs()
        with self._lock:
            job = self._jobs.get(job_id)
            return job.snapshot() if job else None

    def start(self, action: str, params: dict[str, Any] | None = None) -> dict[str, Any]:
        params = dict(params or {})
        if action == "export_hubs" and not params.get("hubs"):
            params.pop("hubs", None)
        needs_pipeline = action in {"build_resource", "build_resources"}
        checker = getattr(
            self.ops,
            "pipeline_available" if needs_pipeline else "hubs_available",
            self.ops.available,
        )
        if self._injected_ops and not checker():
            raise BuildUnavailable(
                "Resource builds need pypath (omnipath-build pipeline)."
                if needs_pipeline
                else "Hub export is not available in this API process."
            )
        stages = self._planned_stages(action, params)
        job = Job(action, params, stages)
        job.log_path = self._job_log_path(job.id)
        self._write_log_header(job)
        try:
            self.job_store.enqueue(job)
        except Exception:
            # The log belongs only to this not-yet-accepted submission.
            job.log_path.unlink(missing_ok=True)
            raise
        self._refresh_jobs()
        return job.snapshot()

    def cancel(self, job_id: str) -> dict[str, Any]:
        return self.job_store.cancel(job_id).snapshot()

    def events(self, job_id: str) -> Iterator[dict[str, Any]]:
        last = None
        while True:
            snapshot = self.get_job(job_id)
            if snapshot is None:
                return
            marker = (snapshot["seq"], snapshot["status"], snapshot.get("cancel_requested"))
            if marker != last:
                yield snapshot
                last = marker
            if snapshot["status"] in TERMINAL_STATUSES:
                return
            time.sleep(0.4)

    def _refresh_jobs(self):
        jobs = self.job_store.list()
        with self._lock:
            self._jobs.update({job.id: job for job in jobs})
            self._current_id = next(
                (job.id for job in jobs if job.status not in TERMINAL_STATUSES),
                jobs[0].id if jobs else None,
            )

    def _recent_jobs(self, limit: int = 50) -> list[Job]:
        jobs = sorted(self._jobs.values(), key=lambda item: item.created_at, reverse=True)
        return jobs[:limit]

    def _job_meta_path(self, job_id: str) -> Path:
        path = self.data_root / "logs" / "jobs" / f"{job_id}.json"
        path.parent.mkdir(parents=True, exist_ok=True)
        return path

    def _persist_job_meta(self, job: Job) -> None:
        self.job_store.save(job)

    def _load_persisted_jobs(self) -> None:
        jobs_dir = self.data_root / "logs" / "jobs"
        if not jobs_dir.is_dir():
            return
        # 1. Load from .json files
        for path in sorted(jobs_dir.glob("*.json"), key=lambda p: p.stat().st_mtime, reverse=True)[
            :50
        ]:
            try:
                data = json.loads(path.read_text(encoding="utf-8"))
                if isinstance(data, dict) and "id" in data:
                    job = Job.from_dict(data)
                    if job.log_path is None:
                        log_file = jobs_dir / f"{job.id}.log"
                        if log_file.is_file():
                            job.log_path = log_file
                    if job.id not in self._jobs:
                        self._jobs[job.id] = job
            except (OSError, ValueError, TypeError):
                logger.exception("Cannot recover job metadata %s", path)
                continue

        # 2. For any .log files that do not have a .json file yet
        for path in sorted(jobs_dir.glob("*.log"), key=lambda p: p.stat().st_mtime, reverse=True)[
            :50
        ]:
            job_id = path.stem
            if job_id in self._jobs:
                continue
            try:
                job = self._parse_job_from_log(path)
                if job is not None and job.id not in self._jobs:
                    self._jobs[job.id] = job
            except (OSError, ValueError, TypeError):
                logger.exception("Cannot recover job metadata %s", path)
                continue

    def _parse_job_from_log(self, path: Path) -> Job | None:
        try:
            text = path.read_text(encoding="utf-8", errors="replace")
            lines = text.splitlines()
            if not lines:
                return None
            action = "build_resources"
            params: dict[str, Any] = {}
            job_id = path.stem
            created_at = path.stat().st_mtime
            for line in lines[:15]:
                if line.startswith("# omnipath admin job"):
                    parts = line.split()
                    if parts:
                        job_id = parts[-1]
                elif line.startswith("# action="):
                    action = line.split("=", 1)[1].strip()
                elif line.startswith("# params="):
                    try:
                        params = json.loads(line.split("=", 1)[1].strip())
                    except (ValueError, TypeError):
                        logger.warning("Malformed legacy log header in %s", path, exc_info=True)
                elif line.startswith("# created="):
                    date_str = line.split("=", 1)[1].strip()
                    try:
                        created_at = time.mktime(time.strptime(date_str, "%Y-%m-%d %H:%M:%S"))
                    except (ValueError, TypeError):
                        logger.warning("Malformed legacy log header in %s", path, exc_info=True)

            stages_seen: list[str] = []
            status = "failed"
            error = "Legacy job log has no completion record"
            for line in lines:
                if "Job finished" in line:
                    status = "done"
                    error = None
                elif "Failed:" in line:
                    status = "failed"
                    error = line.split("Failed:", 1)[1].strip()
                elif "Cancelled" in line:
                    status = "cancelled"
                parts = line.split(": ", 1)
                if len(parts) == 2:
                    first_part = parts[0].strip()
                    subparts = first_part.split()
                    if len(subparts) >= 2:
                        stage_candidate = subparts[-1]
                        if stage_candidate not in stages_seen and not stage_candidate.startswith(
                            "#"
                        ):
                            stages_seen.append(stage_candidate)

            job = Job(action, params, stages_seen)
            job.id = job_id
            job.status = status
            job.error = error
            job.created_at = created_at
            job.started_at = created_at
            job.finished_at = path.stat().st_mtime
            job.logs = lines[-80:]
            job.log_path = path
            for stage in job.stages:
                stage["status"] = "done" if status == "done" else status
            return job
        except (OSError, ValueError, TypeError):
            logger.exception("Cannot recover legacy job log %s", path)
            return None

    def _built_resources(self) -> list[dict[str, Any]]:
        resources_dir = self.data_root / "resources"
        out: list[dict[str, Any]] = []
        if not resources_dir.is_dir():
            return out
        for res_dir in sorted(resources_dir.iterdir()):
            if not res_dir.is_dir():
                continue
            for ver_dir in sorted(res_dir.iterdir()):
                if not ver_dir.is_dir():
                    continue
                entities = _parquet_info(ver_dir / "entity.parquet", "entities")
                relations = _parquet_info(ver_dir / "relation.parquet", "relations")
                payloads = _parquet_info(ver_dir / "evidence_payloads.parquet", "evidence_payloads")
                if not entities["exists"] and not relations["exists"]:
                    continue
                resolution_stats: dict[str, Any] | None = None
                resolution_path = ver_dir / "resolution_stats.json"
                if resolution_path.is_file():
                    try:
                        loaded = json.loads(resolution_path.read_text(encoding="utf-8"))
                        if isinstance(loaded, dict):
                            resolution_stats = loaded
                    except (OSError, ValueError):
                        logger.exception("Cannot read resolution statistics %s", resolution_path)
                        resolution_stats = None
                out.append(
                    {
                        "resource": res_dir.name,
                        "version": ver_dir.name,
                        "entities": entities,
                        "relations": relations,
                        "payloads": payloads,
                        "resolution_stats": resolution_stats,
                        "log": _log_info(ver_dir / "build.log"),
                    }
                )
        return out

    def _available_sources(self) -> list[dict[str, Any]]:
        capabilities = self.job_store.capabilities()
        if capabilities is not None:
            return capabilities.get("sources", [])
        if not self.ops.available():
            return []
        now = time.time()
        cached = self._sources_cache
        if cached is not None and now - cached[0] < 60:
            return cached[1]
        try:
            sources = self.ops.list_sources()
        except (BuildUnavailable, ModuleNotFoundError):
            logger.warning("Optional build source discovery unavailable", exc_info=True)
            return cached[1] if cached is not None else []
        self._sources_cache = (now, sources)
        return sources

    def _planned_stages(self, action: str, params: dict[str, Any]) -> list[str]:
        if action == "export_hubs":
            selected = params.get("hubs")
            names = list(selected) if selected else list(hub_names())
            stages = [f"emit:{name}" for name in names]
            if params.get("include_dictionaries", True):
                stages.append("dictionaries")
            stages.append("manifest")
            return stages
        if action in {"build_resource", "build_resources"}:
            sources = params.get("sources") or ([params["source"]] if params.get("source") else [])
            stages = ["discover"]
            stages.extend(f"build:{src}" for src in sources)
            stages.append("finalize")
            return stages
        raise ValueError(f"Unknown action '{action}'")

    def run_job(self, job: Job, *, should_stop=lambda: False) -> None:
        with self._lock:
            job.status = "running"
            job.started_at = time.time()
            job.seq += 1
            self._cv.notify_all()

        self._persist_job_meta(job)

        def on_progress(event: dict[str, Any]) -> None:
            self._apply_progress(job, event)

        def should_cancel() -> bool:
            persisted = self.job_store.get(job.id)
            job.cancel_requested = (
                job.cancel_requested
                or should_stop()
                or bool(persisted and persisted.cancel_requested)
            )
            return job.cancel_requested

        try:
            result = self._execute(job, on_progress, should_cancel)
            with self._lock:
                now = time.time()
                # Operations honour cancellation by raising BuildCancelled. A normal
                # return means the work (e.g. a resource publication) completed.
                if job.cancel_requested:
                    self._append_log(job, "Cancel requested after the work had completed")
                job.status = "done"
                job.result = _jsonable(result)
                self._close_open_stages(job, "done", "", now)
                self._append_log(job, "Job finished")
                job.finished_at = now
            self._persist_job_logs(job)
            self._persist_job_meta(job)
            with self._lock:
                job.seq += 1
                self._cv.notify_all()
        except Exception as exc:
            cancelled = exc.__class__.__name__ == "BuildCancelled" or job.cancel_requested
            with self._lock:
                now = time.time()
                job.finished_at = now
                if cancelled:
                    job.status = "cancelled"
                    self._close_open_stages(job, "cancelled", "Cancelled", now)
                    self._append_log(job, "Cancelled")
                else:
                    job.status = "failed"
                    job.error = str(exc)
                    if isinstance(getattr(exc, "results", None), dict):
                        job.result = _jsonable(exc.results)
                    self._close_open_stages(job, "failed", str(exc), now)
                    self._append_log(job, f"Failed: {exc}")
                    self._append_log(job, traceback.format_exc()[-1500:])
            self._persist_job_logs(job)
            self._persist_job_meta(job)
            with self._lock:
                job.seq += 1
                self._cv.notify_all()

    def _execute(
        self,
        job: Job,
        on_progress: Callable[[dict[str, Any]], None],
        should_cancel: Callable[[], bool],
    ) -> Any:
        hubs_dir = self.hubs_dir()
        if job.action == "export_hubs":
            max_records = job.params.get("max_records")
            if max_records == 0:
                max_records = None
            parallel = int(job.params.get("parallel") or 1)
            return self.ops.export_hubs(
                output_dir=hubs_dir,
                hubs=job.params.get("hubs"),
                max_records=max_records,
                include_dictionaries=job.params.get("include_dictionaries", True),
                parallel=parallel,
                on_progress=on_progress,
                should_cancel=should_cancel,
            )
        if job.action in {"build_resource", "build_resources"}:
            sources = job.params.get("sources") or [job.params["source"]]
            max_records = job.params.get("max_records")
            if max_records == 0:
                max_records = None
            parallel = int(job.params.get("parallel") or 1)
            return self.ops.build_all(
                sources=sources,
                version=job.params.get("version"),
                versions=job.params.get("versions"),
                output_dir=self.data_root,
                max_records=max_records,
                batch_size=int(job.params.get("batch_size") or 50000),
                force_refresh=bool(job.params.get("force_refresh") or False),
                progress=False,
                parallel=parallel,
                on_progress=on_progress,
                should_cancel=should_cancel,
            )
        raise ValueError(f"Unknown action '{job.action}'")

    def _apply_progress(self, job: Job, event: dict[str, Any]) -> None:
        stage_id = str(event.get("stage") or "unknown")
        status = str(event.get("status") or "running")
        message = str(event.get("message") or "")
        with self._lock:
            now = time.time()
            stage = next((item for item in job.stages if item["id"] == stage_id), None)
            if stage is None:
                stage = _make_stage(stage_id)
                finalize = next(
                    (index for index, item in enumerate(job.stages) if item["id"] == "finalize"),
                    None,
                )
                if finalize is None:
                    job.stages.append(stage)
                else:
                    job.stages.insert(finalize, stage)
            _set_stage_status(stage, status, now)
            if message:
                stage["message"] = message
            if "current" in event:
                stage["current"] = event.get("current")
            if "total" in event:
                stage["total"] = event.get("total")
            if event.get("counts"):
                stage["counts"] = dict(event["counts"])
            if message:
                self._append_log(job, f"{stage_id}: {message}")
            job.seq += 1
            self._cv.notify_all()
            self._persist_job_meta(job)

    def _close_open_stages(self, job: Job, status: str, message: str, now: float) -> None:
        for stage in job.stages:
            if stage["status"] == "running":
                _set_stage_status(stage, status, now)
                if message:
                    stage["message"] = message

    def _append_log(self, job: Job, line: str) -> None:
        stamp = time.strftime("%H:%M:%S")
        text = f"{stamp} {line}"
        job.logs.append(text)
        if len(job.logs) > 2000:
            job.logs = job.logs[-2000:]
        self._write_log_line(job, text)

    def _job_log_path(self, job_id: str) -> Path:
        path = self.data_root / "logs" / "jobs" / f"{job_id}.log"
        path.parent.mkdir(parents=True, exist_ok=True)
        return path

    def _hub_log_path(self, name: str) -> Path:
        return self.hubs_dir() / "logs" / f"{name}.log"

    def _resource_alias_log_path(self, name: str) -> Path:
        return self.data_root / "logs" / "resources" / f"{name}.log"

    def _resource_log_map(self) -> dict[str, dict[str, Any]]:
        logs: dict[str, dict[str, Any]] = {}
        alias_dir = self.data_root / "logs" / "resources"
        if alias_dir.is_dir():
            for path in alias_dir.glob("*.log"):
                logs[path.stem] = _log_info(path)
        return logs

    def _write_log_header(self, job: Job) -> None:
        if job.log_path is None:
            return
        lines = [
            f"# omnipath admin job {job.id}",
            f"# action={job.action}",
            f"# params={json.dumps(job.params, default=str, sort_keys=True)}",
            f"# created={time.strftime('%Y-%m-%d %H:%M:%S')}",
            "",
        ]
        job.log_path.parent.mkdir(parents=True, exist_ok=True)
        job.log_path.write_text("\n".join(lines), encoding="utf-8")

    def _write_log_line(self, job: Job, text: str) -> None:
        if job.log_path is None:
            return
        try:
            job.log_path.parent.mkdir(parents=True, exist_ok=True)
            with job.log_path.open("a", encoding="utf-8") as handle:
                handle.write(text + "\n")
        except OSError:
            logger.exception("Cannot append job log %s", job.log_path)
            raise

    def _persist_job_logs(self, job: Job) -> None:
        source = job.log_path
        if source is None or not source.is_file():
            return
        for target in self._log_targets(job):
            try:
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(source, target)
            except OSError:
                logger.exception("Cannot persist build log %s", target)
                raise

    def _log_targets(self, job: Job) -> list[Path]:
        paths: list[Path] = []
        if job.action == "export_hubs":
            names = list(job.params.get("hubs") or hub_names())
            paths.extend(self._hub_log_path(name) for name in names)
        elif job.action in {"build_resource", "build_resources"}:
            sources = job.params.get("sources") or (
                [job.params["source"]] if job.params.get("source") else []
            )
            paths.extend(self._resource_alias_log_path(str(name)) for name in sources)
            result = job.result if isinstance(job.result, dict) else {}
            for key, info in result.items():
                resource = ""
                version = ""
                if isinstance(info, dict):
                    resource = str(info.get("resource") or "")
                    version = str(info.get("version") or "")
                if not resource and isinstance(key, str) and "/" in key:
                    resource, version = key.split("/", 1)
                if resource and version:
                    paths.append(self.data_root / "resources" / resource / version / "build.log")
        return paths

    def read_log(self, kind: str, name: str, version: str | None = None) -> dict[str, Any]:
        if not _SAFE_LOG_NAME.match(kind) or not _SAFE_LOG_NAME.match(name):
            raise ValueError("Invalid log name")
        if version is not None and not _SAFE_LOG_NAME.match(version):
            raise ValueError("Invalid log version")
        path = self._resolve_log_path(kind, name, version)
        payload = {
            "kind": kind,
            "name": name,
            "version": version,
            "path": str(path) if path else None,
            "exists": bool(path and path.is_file()),
            "mtime": path.stat().st_mtime if path and path.is_file() else None,
            "text": "",
        }
        if path is not None and path.is_file():
            payload["text"] = _read_log_text(path)
        return payload

    def _resolve_log_path(self, kind: str, name: str, version: str | None) -> Path | None:
        if kind == "job":
            with self._lock:
                job = self._jobs.get(name)
                if job is not None and job.log_path is not None:
                    return job.log_path
            candidate = self.data_root / "logs" / "jobs" / f"{name}.log"
            return candidate if candidate.is_file() else None
        if kind == "hub":
            return self._hub_log_path(name)
        if kind == "resource":
            if version:
                return self.data_root / "resources" / name / version / "build.log"
            alias = self._resource_alias_log_path(name)
            if alias.is_file():
                return alias
            built = [item for item in self._built_resources() if item["resource"] == name]
            usable = [item for item in built if item.get("log", {}).get("exists")]
            if usable:
                latest = max(usable, key=lambda item: item["log"].get("mtime") or 0)
                return Path(latest["log"]["path"])
            return alias
        raise ValueError(f"Unknown log kind '{kind}'")

    def delete_resource_version(self, resource: str, version: str) -> dict[str, Any]:
        if not _SAFE_LOG_NAME.match(resource) or not _SAFE_LOG_NAME.match(version):
            raise ValueError("Invalid resource or version name")
        references = ReleaseStore(self.data_root).references(resource, version)
        if references:
            raise ValueError(
                f"Resource {resource}/{version} is pinned by OmniPath {', '.join(references)}"
            )
        target_dir = self.data_root / "resources" / resource / version
        if not target_dir.is_dir():
            raise FileNotFoundError(f"Version {version} of resource {resource} not found")
        shutil.rmtree(target_dir)
        res_dir = self.data_root / "resources" / resource
        if res_dir.is_dir() and not any(res_dir.iterdir()):
            try:
                res_dir.rmdir()
            except OSError:
                logger.warning(
                    "Could not remove empty resource directory %s", res_dir, exc_info=True
                )
        return {"status": "deleted", "resource": resource, "version": version}


def _read_log_text(path: Path) -> str:
    return read_log_tail(path)


def _jsonable(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    return value


def register_admin_routes(app):
    """Compatibility import for consumers constructing their own app."""
    from omnipath_api.routers.admin import register_admin_routes as register

    register(app)
