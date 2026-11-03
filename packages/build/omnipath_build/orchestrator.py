"""Memory-aware admission of isolated resource builds under a shared CPU budget.

Linux/systemd gives the workers an aggregate cgroup ceiling and an individual
memory ceiling per attempt. CPUs are shared work-conservingly by the slice.
The local executor is explicitly advisory, for development and portable tests.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from datetime import datetime, timezone
import json
import fcntl
import math
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import sys
import time
import uuid
from typing import Any

import psutil

from omnipath_core.versioning import RESOURCE_FILES, SERVING_FILES

from .versioning import validate_source, validate_version

GIB = 1024**3
MIB = 1024**2


def memory_bytes(value: str) -> int:
    if not isinstance(value, str):
        raise ValueError("RAM must be a size string such as 10GiB")
    match = re.fullmatch(r"\s*(\d+(?:\.\d*)?|\.\d+)\s*([KMGT]?)(I?B)\s*", value.upper())
    if not match:
        raise ValueError("RAM must be a positive size such as 10GiB or 8000MB")
    number, prefix, unit = match.groups()
    size = int(
        float(number)
        * (1024 if unit == "IB" else 1000) ** ("KMGT".index(prefix) + 1 if prefix else 0)
    )
    if size <= 0:
        raise ValueError("RAM must be positive")
    return size


def atomic_json(path: Path, data: Any):
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(data, default=str, indent=2) + "\n")
    temp.replace(path)


def update_history(path: Path, key: str, value: dict):
    """Concurrent scheduler runs merge learning without losing each other's entries."""
    with path.with_suffix(".lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        try:
            current = json.loads(path.read_text()) if path.exists() else {}
            current[key] = value
            atomic_json(path, current)
        finally:
            fcntl.flock(lock, fcntl.LOCK_UN)


@dataclass(frozen=True)
class Budget:
    ram: int
    cpus: int
    jobs: int | None = None
    worker_ram: int = 4 * GIB
    reserve_cpus: int = 1

    def __post_init__(self):
        if self.ram < 256 * MIB or self.worker_ram < 128 * MIB:
            raise ValueError("Use at least 256MiB total RAM and 128MiB per worker")
        if self.reserve_cpus < 1:
            raise ValueError("At least one CPU must remain outside the build")
        if self.cpus < 1 or (self.jobs is not None and self.jobs < 1):
            raise ValueError("CPUs and jobs must be positive integers")

    @property
    def max_jobs(self):
        return min(self.jobs or self.cpus, self.cpus)


def cpu_allocation(budget: Budget):
    available = (
        sorted(os.sched_getaffinity(0))
        if hasattr(os, "sched_getaffinity")
        else list(range(os.cpu_count() or 1))
    )
    count = min(budget.cpus, len(available) - budget.reserve_cpus)
    if count < 1:
        raise ValueError("Not enough CPUs to build while preserving the requested CPU headroom")
    return replace(budget, cpus=count), available[:count]


@dataclass
class Job:
    source: str
    version: str
    ram: int
    duration: float = 0
    attempts: int = 0
    status: str = "queued"
    message: str = "Waiting for memory"
    peak: int = 0
    fraction: float | None = None
    result: dict = field(default_factory=dict)
    workers: int | None = None


@dataclass
class Running:
    job: Job
    process: subprocess.Popen
    directory: Path
    log: Any
    events: Any
    unit: str = ""
    started: float = field(default_factory=time.monotonic)
    ram: int = 0
    cpu: float = 0
    memory_failure: bool = False


class LocalExecutor:
    """Advisory reservations, without an OS-level total RAM or CPU ceiling."""

    name = "local (advisory limits)"

    def open(self, run_id: str, budget: Budget):
        pass

    def start(self, command, log, job, attempt):
        return subprocess.Popen(
            command, stdout=log, stderr=subprocess.STDOUT, start_new_session=True
        ), ""

    def sample(self, running):
        try:
            root = psutil.Process(running.process.pid)
            try:
                processes = [root, *root.children(recursive=True)]
            except (psutil.AccessDenied, PermissionError):
                # Process enumeration may be restricted in development sandboxes.
                # Local reservations are advisory; retain parent sampling when possible.
                processes = [root]
                if not getattr(self, "_warned_process_access", False):
                    import logging

                    logging.getLogger(__name__).warning(
                        "Cannot enumerate worker descendants; local resource metrics cover the parent only"
                    )
                    self._warned_process_access = True
            memory = cpu = 0
            for process in processes:
                try:
                    memory += process.memory_info().rss
                    times = process.cpu_times()
                    cpu += times.user + times.system
                except (psutil.NoSuchProcess, psutil.AccessDenied, PermissionError):
                    pass
            running.ram, running.cpu = memory, cpu
            running.job.peak = max(running.job.peak, memory)
        except psutil.NoSuchProcess:
            pass

    def stop(self, running, force=False):
        try:
            os.killpg(running.process.pid, signal.SIGKILL if force else signal.SIGTERM)
        except ProcessLookupError:
            pass

    def finish(self, running):
        self.stop(running, force=True)

    def close(self):
        pass


class SystemdExecutor(LocalExecutor):
    """Unique runtime slice, with independently killable per-attempt scopes."""

    name = "systemd (enforced worker limits)"

    def __init__(self, user=False):
        self.user = ["--user"] if user else []
        self.slice = ""

    def ctl(self, *args, check=True):
        return subprocess.run(
            ["systemctl", *self.user, *args],
            check=check,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=20,
        )

    def open(self, run_id, budget):
        if (
            not sys.platform.startswith("linux")
            or not Path("/sys/fs/cgroup/cgroup.controllers").exists()
        ):
            raise ValueError(
                "Enforced scheduling requires Linux cgroup v2 and systemd; use --executor local only for advisory development runs"
            )
        self.slice = f"omnipathbuild{run_id}.slice"
        try:
            self.ctl("start", self.slice)
            self.ctl(
                "set-property",
                "--runtime",
                self.slice,
                f"MemoryMax={budget.ram}",
                "MemorySwapMax=0",
                f"CPUQuota={budget.cpus * 100}%",
                "MemoryAccounting=yes",
                "CPUAccounting=yes",
                "AllowedCPUs=" + ",".join(str(cpu) for cpu in self.allowed_cpus),
            )
        except Exception:
            self.close()
            raise

    def start(self, command, log, job, attempt):
        unit = f"omnipathbuild{uuid.uuid4().hex}.scope"
        args = [
            "systemd-run",
            *self.user,
            "--quiet",
            "--scope",
            f"--unit={unit}",
            f"--slice={self.slice}",
            f"--property=MemoryMax={job.ram}",
            "--property=MemorySwapMax=0",
            "--property=OOMPolicy=kill",
            "--property=CPUWeight=100",
            "--property=MemoryAccounting=yes",
            "--property=CPUAccounting=yes",
            *command,
        ]
        return subprocess.Popen(
            args, stdout=log, stderr=subprocess.STDOUT, start_new_session=True
        ), unit

    def sample(self, running):
        props = self.ctl(
            "show",
            running.unit,
            "-p",
            "MemoryCurrent",
            "-p",
            "MemoryPeak",
            "-p",
            "CPUUsageNSec",
            "-p",
            "Result",
            check=False,
        )
        values = dict(line.split("=", 1) for line in props.stdout.splitlines() if "=" in line)

        def number(key):
            value = values.get(key, "")
            return int(value) if value.isdecimal() and int(value) < 2**63 else 0

        running.ram = number("MemoryCurrent")
        running.cpu = max(running.cpu, number("CPUUsageNSec") / 1e9)
        running.job.peak = max(running.job.peak, running.ram, number("MemoryPeak"))
        running.memory_failure |= values.get("Result") == "oom-kill"

    def stop(self, running, force=False):
        self.ctl(
            "kill",
            "--kill-whom=all",
            f"--signal={'SIGKILL' if force else 'SIGTERM'}",
            running.unit,
            check=False,
        )
        # Also reap systemd-run, including a launcher that failed before scope creation.
        if force:
            super().stop(running, force=True)

    def finish(self, running):
        # No subprocess may outlive its reservation, even after its Python parent exits.
        self.ctl("stop", running.unit, check=False)
        self.ctl("reset-failed", running.unit, check=False)

    def close(self):
        if self.slice:
            self.ctl("stop", self.slice, check=False)
            self.ctl("reset-failed", self.slice, check=False)
            self.ctl("revert", self.slice, check=False)


class ConsoleProgress:
    """TTY dashboard; throttled, timestamped lines when redirected."""

    def __init__(self, stream=None):
        self.stream = stream or sys.stderr
        self.tty = self.stream.isatty()
        self.lines = 0
        self.last = {}
        self.last_time = 0

    def __call__(self, snapshot):
        now = time.monotonic()
        jobs = snapshot["jobs"]
        elapsed = int(snapshot["elapsed_seconds"])
        counts = {
            state: sum(j["status"] == state for j in jobs)
            for state in ("queued", "running", "success", "failed", "cancelled", "skipped")
        }
        title = (
            f"OmniPath  {elapsed // 60:02d}:{elapsed % 60:02d}  |  "
            f"{counts['success']} done  {counts['running']} running  {counts['queued']} queued  {counts['skipped']} skipped  {counts['failed']} failed"
        )
        usage = (
            f"RAM {snapshot['ram_bytes'] / GIB:.2f} GiB used / {snapshot['reserved_bytes'] / GIB:.2f} reserved / "
            f"{snapshot['budget']['ram'] / GIB:.2f} budget   CPU {snapshot['cpu_cores']:.1f} / {snapshot['budget']['cpus']}"
        )

        def clean(value):
            return re.sub(r"[\x00-\x1f\x7f-\x9f]", " ", str(value))

        if self.tty:
            width = max(40, shutil.get_terminal_size().columns - 1)
            lines = [title, usage, snapshot["executor"]]
            if snapshot.get("error"):
                lines.append(snapshot["error"])
            for job in jobs:
                if job["status"] == "running":
                    fraction = job.get("fraction")
                    bar = (
                        "··········"
                        if fraction is None
                        else "━" * int(fraction * 10) + "─" * (10 - int(fraction * 10))
                    )
                    lines.append(
                        f"  ▶ {job['source']:<20} {job['ram_bytes'] / GIB:5.2f} GiB  [{bar}] {job['message']}"
                    )
            if not counts["running"]:
                lines.append(
                    "  ✓ Finished" if snapshot["status"] == "success" else f"  {snapshot['status']}"
                )
            prefix = f"\x1b[{self.lines}F" if self.lines else ""
            self.stream.write(
                prefix + "\x1b[J" + "\n".join(clean(line)[:width] for line in lines) + "\n"
            )
            self.lines = len(lines)
        else:
            if not self.last_time or now - self.last_time >= 30 or snapshot["status"] != "running":
                print(title + " | " + usage + " | " + snapshot["executor"], file=self.stream)
                self.last_time = now
                if snapshot.get("error"):
                    print(clean(snapshot["error"]), file=self.stream)
            for job in jobs:
                previous, timestamp = self.last.get(job["source"], (None, 0))
                state = (job["status"], job["message"])
                if state != previous and (
                    previous is None or state[0] != previous[0] or now - timestamp >= 10
                ):
                    print(
                        f"{elapsed // 60:02d}:{elapsed % 60:02d}  {job['status'].upper():9} {job['source']:<20} {clean(job['message'])}",
                        file=self.stream,
                    )
                    self.last[job["source"]] = state, now
        self.stream.flush()


def fit_jobs(pending: list[Job], running: list[Running], budget: Budget) -> list[Job]:
    """Longest known jobs first; backfill smaller jobs rather than block the queue."""
    free = budget.ram - sum(item.job.ram for item in running)
    slots = budget.max_jobs - len(running)
    selected = []
    for job in sorted(pending, key=lambda job: (-job.duration, -job.ram, job.source)):
        if slots and job.ram <= free:
            selected.append(job)
            free -= job.ram
            slots -= 1
    return selected


def _read_events(item):
    while True:
        position = item.events.tell()
        line = item.events.readline()
        if not line:
            break
        if not line.endswith("\n"):
            item.events.seek(position)
            break
        try:
            event = json.loads(line)
            item.job.message = str(event.get("message") or event.get("stage") or "Building")
            total = event.get("total")
            item.job.fraction = (
                min(1.0, max(0.0, float(event.get("current", 0)) / total))
                if isinstance(total, (int, float)) and total > 0
                else None
            )
        except (ValueError, TypeError):
            pass


def orchestrate(
    sources,
    *,
    budget: Budget,
    version=None,
    versions=None,
    output_dir="data",
    max_records=None,
    batch_size=50000,
    executor=None,
    on_progress=None,
    should_cancel=None,
    skip_existing=False,
    profiles=None,
    max_retries=2,
    poll_interval=1.0,
    library_dir=None,
    worker_command=None,
    min_free_disk=0,
):
    """Run selected immutable resource builds. Return durable summary, including failures.

    `worker_command` is a test seam implementing the spec/events/result protocol.
    It is deliberately not exposed as a shell command option in the public CLI.
    """
    if min_free_disk < 0:
        raise ValueError("Free disk reserve cannot be negative")
    if (
        max_retries < 0
        or poll_interval <= 0
        or batch_size < 1
        or (max_records is not None and max_records < 1)
    ):
        raise ValueError("Invalid retry, polling, batch or record limit")
    if versions is not None and not isinstance(versions, dict):
        raise ValueError("Versions must be a JSON object mapping resources to numeric versions")
    if profiles is not None and (
        not isinstance(profiles, dict) or any(not isinstance(p, dict) for p in profiles.values())
    ):
        raise ValueError("Profiles must map resources to objects containing a ram size")
    requested_cpus = budget.cpus
    budget, allowed_cpus = cpu_allocation(budget)
    root = Path(output_dir).resolve()
    names = list(dict.fromkeys(validate_source(source.strip().lower()) for source in sources))
    if not names:
        raise ValueError("Select at least one resource")
    selected = {name: validate_version((versions or {}).get(name, version)) for name in names}
    history_path = root / ".build-history.json"
    try:
        history = json.loads(history_path.read_text())
    except FileNotFoundError:
        history = {}

    def profile_key(source):
        return f"{source}:{max_records if max_records is not None else 'all'}"

    jobs = []
    for name in names:
        target = root / "resources" / name / selected[name]
        if target.exists():
            if not skip_existing:
                raise FileExistsError(
                    f"{target} already exists; use a new version or --skip-existing"
                )
            manifest = json.loads((target / "build_manifest.json").read_text())
            if (
                manifest.get("resource") != name
                or manifest.get("version") != selected[name]
                or not all((target / file).is_file() for file in RESOURCE_FILES + SERVING_FILES)
            ):
                raise ValueError(f"Cannot skip incomplete resource version {target}")
        estimate = history.get(profile_key(name), {})
        explicit = (profiles or {}).get(name, {})
        ram = (
            memory_bytes(explicit["ram"])
            if "ram" in explicit
            else min(
                budget.ram,
                max(128 * MIB, int(estimate.get("reservation_bytes", budget.worker_ram))),
            )
        )
        if ram > budget.ram or ram < 128 * MIB:
            raise ValueError(
                f"{name}: memory request must be between 128MiB and the total RAM budget"
            )
        jobs.append(
            Job(
                name,
                selected[name],
                ram,
                float(estimate.get("elapsed_seconds", 0)),
                status="skipped" if target.exists() else "queued",
                message="Existing complete version" if target.exists() else "Waiting for memory",
            )
        )
    run_id = uuid.uuid4().hex
    run_dir = (
        root / ".build-runs" / (datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ-") + run_id[:8])
    )
    run_dir.mkdir(parents=True)
    executor = executor or SystemdExecutor()
    executor.allowed_cpus = allowed_cpus
    pending = [job for job in jobs if job.status == "queued"]
    running = []
    locks = {}
    started = time.monotonic()
    cpu_previous = (started, 0.0)
    finished_cpu = 0.0
    status = "running"
    summary = {}
    cancelled = False
    run_error = None

    def publish_snapshot():
        nonlocal cpu_previous, summary
        now = time.monotonic()
        cpu = finished_cpu + sum(item.cpu for item in running)
        cores = max(0.0, cpu - cpu_previous[1]) / max(0.001, now - cpu_previous[0])
        cpu_previous = now, cpu
        summary = {
            "run_id": run_id,
            "run_dir": str(run_dir),
            "status": status,
            "executor": executor.name,
            "budget": {
                "ram": budget.ram,
                "cpus": budget.cpus,
                "requested_cpus": requested_cpus,
                "reserve_cpus": budget.reserve_cpus,
                "allowed_cpus": allowed_cpus,
            },
            "elapsed_seconds": now - started,
            "cpu_cores": cores,
            "error": run_error,
            "min_free_disk_bytes": min_free_disk,
            "ram_bytes": sum(item.ram for item in running),
            "reserved_bytes": sum(item.job.ram for item in running),
            "jobs": [
                {
                    "source": job.source,
                    "version": job.version,
                    "status": job.status,
                    "attempts": job.attempts,
                    "reservation_bytes": job.ram,
                    "peak_bytes": job.peak,
                    "ram_bytes": next((r.ram for r in running if r.job is job), 0),
                    "message": job.message,
                    "fraction": job.fraction,
                    "result": job.result,
                }
                for job in jobs
            ],
        }
        atomic_json(run_dir / "summary.json", summary)
        if on_progress:
            on_progress(summary)

    try:
        executor.open(run_id, budget)
        while pending or running:
            if min_free_disk and shutil.disk_usage(root).free < min_free_disk:
                run_error = f"Free disk space fell below the {min_free_disk / GIB:.1f} GiB reserve; stopping workers"
                break
            if should_cancel and should_cancel():
                cancelled = True
                break
            for job in fit_jobs(pending, running, budget):
                pending.remove(job)
                lock = root / "resources" / job.source / f".{job.version}.lock"
                if job.source not in locks:
                    lock.parent.mkdir(parents=True, exist_ok=True)
                    try:
                        from .locking import BuildLock

                        locks[job.source] = BuildLock(lock).__enter__()
                        if (lock.parent / job.version).exists():
                            raise FileExistsError("Resource version was published while queued")
                    except FileExistsError as exc:
                        job.status, job.message = "failed", str(exc)
                        continue
                job.attempts += 1
                directory = run_dir / job.source / f"attempt-{job.attempts}"
                directory.mkdir(parents=True)
                # Divide capacity when resource concurrency is explicitly requested.
                # Preparation workers use the CPU budget, then exit before finalization.
                resource_cpus = max(1, budget.cpus // budget.max_jobs)
                if job.workers is None:
                    job.workers = max(1, min(resource_cpus, job.ram // GIB))
                batch_workers = job.workers
                batch_share = job.ram // batch_workers
                spec = {
                    "cpus": resource_cpus,
                    "allowed_cpus": allowed_cpus,
                    "tmpdir": str(directory / "tmp"),
                    "events": str(directory / "events.jsonl"),
                    "result": str(directory / "result.json"),
                    "build": {
                        "source": job.source,
                        "version": job.version,
                        "output_dir": str(directory / "output"),
                        "max_records": max_records,
                        "batch_size": batch_size,
                        "max_batch_records": batch_size,
                        "batch_workers": batch_workers,
                        "resource_ram_bytes": job.ram,
                        "min_free_disk_bytes": min_free_disk,
                        "raw_batch_bytes": batch_share,
                        "max_batch_bytes": min(batch_share, 64 * MIB),
                        "library_dir": str(Path(library_dir).resolve())
                        if library_dir
                        else os.environ.get(
                            "OMNIPATH_LIBRARY_DIR", str(root / "reference/library")
                        ),
                    },
                }
                atomic_json(directory / "spec.json", spec)
                Path(spec["events"]).touch()
                log = (directory / "console.log").open("w")
                events = Path(spec["events"]).open()
                command = [
                    *(
                        worker_command
                        or [sys.executable, "-m", "omnipath_build.orchestrator_worker"]
                    ),
                    str(directory / "spec.json"),
                ]
                try:
                    process, unit = executor.start(command, log, job, directory)
                except Exception as exc:
                    log.close()
                    events.close()
                    job.status, job.message = "failed", f"Worker could not start: {exc}"
                    continue
                running.append(Running(job, process, directory, log, events, unit))
                job.status, job.message = "running", f"Starting attempt {job.attempts}"
            for item in list(running):
                executor.sample(item)
                _read_events(item)
                if item.process.poll() is None:
                    continue
                executor.finish(item)
                item.log.close()
                item.events.close()
                running.remove(item)
                job = item.job
                result_path = item.directory / "result.json"
                try:
                    result = json.loads(result_path.read_text())
                except (FileNotFoundError, ValueError):
                    result = {
                        "status": "failed",
                        "error": f"Worker exited {item.process.returncode}; see {item.directory / 'console.log'}",
                    }
                job.peak = max(job.peak, int(result.get("peak_rss_bytes", 0)))
                finished_cpu += max(item.cpu, float(result.get("cpu_seconds", 0)))
                success = item.process.returncode == 0 and result.get("status") == "success"
                if success:
                    built = item.directory / "output/resources" / job.source / job.version
                    target = root / "resources" / job.source / job.version
                    try:
                        manifest = json.loads((built / "build_manifest.json").read_text())
                        if (
                            manifest.get("resource") != job.source
                            or manifest.get("version") != job.version
                        ):
                            raise ValueError("Worker manifest identity does not match its job")
                        for name in RESOURCE_FILES + SERVING_FILES:
                            if not (built / name).is_file():
                                raise ValueError(f"Worker output is missing {name}")
                        if target.exists():
                            raise FileExistsError(f"{target} already exists")
                        built.rename(target)
                        result["files"] = {
                            name: str(target / Path(path).name)
                            for name, path in result.get("files", {}).items()
                        }
                        for key in ("resolution_stats_path", "manifest_path"):
                            if key in result:
                                result[key] = str(target / Path(result[key]).name)
                        job.status, job.message = (
                            "success",
                            f"Finished in {result.get('elapsed_seconds', 0):.1f}s",
                        )
                        history[profile_key(job.source)] = {
                            "reservation_bytes": min(
                                budget.ram, max(512 * MIB, math.ceil(job.peak * 1.3))
                            ),
                            "peak_bytes": job.peak,
                            "elapsed_seconds": result.get("elapsed_seconds", 0),
                        }
                        try:
                            update_history(
                                history_path,
                                profile_key(job.source),
                                history[profile_key(job.source)],
                            )
                        except (OSError, ValueError) as exc:
                            result["history_warning"] = str(exc)
                    except Exception as exc:
                        job.status, job.message = "failed", f"Publication failed: {exc}"
                else:
                    oom = item.memory_failure or result.get("memory_failure", False)
                    result["memory_failure"] = bool(oom)
                    if oom:
                        result["error"] = "Memory limit exceeded: " + str(
                            result.get("error", "worker stopped")
                        )
                    # A retry doubles the memory and keeps the workers, so each worker's
                    # share doubles; only at the budget cap are the workers halved instead.
                    next_ram = min(budget.ram, job.ram * 2)
                    workers = job.workers or 1
                    next_workers = workers if next_ram > job.ram else max(1, workers // 2)
                    if (
                        oom
                        and job.attempts <= max_retries
                        and (next_ram > job.ram or next_workers < workers)
                    ):
                        job.ram, job.workers = next_ram, next_workers
                        job.status, job.message = (
                            "queued",
                            f"Memory retry with {job.ram / GIB:.2f} GiB",
                        )
                        pending.append(job)
                    else:
                        job.status, job.message = (
                            "failed",
                            str(result.get("error", "Worker failed")),
                        )
                job.result = result
                # Only this attempt's owned workspace is removed; diagnostics remain.
                shutil.rmtree(item.directory / "output", ignore_errors=True)
                shutil.rmtree(item.directory / "tmp", ignore_errors=True)
            publish_snapshot()
            if pending or running:
                time.sleep(poll_interval)
        status = (
            "cancelled"
            if cancelled
            else ("failed" if run_error or any(j.status == "failed" for j in jobs) else "success")
        )
    except KeyboardInterrupt:
        cancelled, status = True, "cancelled"
    except Exception:
        status = "failed"
        raise
    finally:
        for item in running:
            executor.stop(item)
        deadline = time.monotonic() + 5
        for item in running:
            try:
                item.process.wait(timeout=max(0.01, deadline - time.monotonic()))
            except subprocess.TimeoutExpired:
                executor.stop(item, force=True)
                item.process.wait(timeout=5)
            executor.finish(item)
            item.log.close()
            item.events.close()
            item.job.status, item.job.message = "cancelled", run_error or "Worker stopped"
            shutil.rmtree(item.directory / "output", ignore_errors=True)
            shutil.rmtree(item.directory / "tmp", ignore_errors=True)
        running.clear()
        for job in pending:
            job.status, job.message = "cancelled", "Not started"
        for lock in locks.values():
            lock.close()
        executor.close()
        publish_snapshot()
    return summary
