import io
import json
import os
from pathlib import Path
import sys
import time
from types import SimpleNamespace

import psutil
import pytest

from omnipath_build.locking import BuildLock


from omnipath_build.orchestrator import (
    Budget,
    ConsoleProgress,
    GIB,
    MIB,
    Job,
    LocalExecutor,
    SystemdExecutor,
    cpu_allocation,
    fit_jobs,
    memory_bytes,
    orchestrate,
)

WORKER = [sys.executable, str(Path(__file__).parent / "fixtures/scheduler_worker.py")]


def assert_locks_released(root):
    for path in root.glob("resources/*/.*.lock"):
        with BuildLock(path):
            pass


def test_memory_units_and_required_cpu_headroom(monkeypatch):
    assert memory_bytes("10GiB") == 10 * GIB
    assert memory_bytes("1.5GB") == 1500000000
    for value in ("all", "-1GB", "0B", "nanGB"):
        with pytest.raises(ValueError):
            memory_bytes(value)
    monkeypatch.setattr(os, "sched_getaffinity", lambda pid: {3, 5, 8, 9}, raising=False)
    budget, cpus = cpu_allocation(Budget(GIB, 100, reserve_cpus=1))
    assert budget.cpus == 3 and cpus == [3, 5, 8]
    with pytest.raises(ValueError, match="At least one"):
        Budget(GIB, 2, reserve_cpus=0)
    with pytest.raises(ValueError, match="Not enough"):
        cpu_allocation(Budget(GIB, 2, reserve_cpus=4))


def test_scheduler_backfills_without_overcommitting():
    budget = Budget(10 * GIB, 4)
    active = [SimpleNamespace(job=Job("active", "1", 6 * GIB))]
    pending = [
        Job("large", "1", 8 * GIB, 100),
        Job("small", "1", 3 * GIB, 10),
        Job("tiny", "1", GIB),
    ]
    assert [j.source for j in fit_jobs(pending, active, budget)] == ["small", "tiny"]


def run_fixture(tmp_path, sources, **kwargs):
    return orchestrate(
        sources,
        budget=Budget(GIB, 2, worker_ram=512 * MIB),
        version="1",
        output_dir=tmp_path,
        executor=LocalExecutor(),
        poll_interval=0.02,
        worker_command=WORKER,
        **kwargs,
    )


def test_processes_publish_retry_continue_and_record_history(tmp_path):
    snapshots = []
    result = run_fixture(
        tmp_path, ["first", "second", "retry", "fail"], on_progress=snapshots.append
    )
    assert result["status"] == "failed"
    jobs = {job["source"]: job for job in result["jobs"]}
    assert jobs["retry"]["attempts"] == 2 and jobs["retry"]["status"] == "success"
    assert jobs["fail"]["status"] == "failed"
    assert jobs["first"]["status"] == jobs["second"]["status"] == "success"
    assert max(snapshot["reserved_bytes"] for snapshot in snapshots) <= GIB
    assert (
        max(sum(j["status"] == "running" for j in snapshot["jobs"]) for snapshot in snapshots) == 2
    )
    assert_locks_released(tmp_path)
    assert (tmp_path / "resources/retry/1/build_manifest.json").exists()
    history = json.loads((tmp_path / ".build-history.json").read_text())
    assert "retry:all" in history and "fail:all" not in history
    assert list(Path(result["run_dir"]).glob("retry/attempt-*/console.log"))
    specs = [
        json.loads(p.read_text())
        for p in sorted(Path(result["run_dir"]).glob("retry/attempt-*/spec.json"))
    ]
    assert specs[1]["build"]["resource_ram_bytes"] > specs[0]["build"]["resource_ram_bytes"]
    assert not list(Path(result["run_dir"]).glob("*/attempt-*/output"))


def test_resume_and_existing_version_preflight(tmp_path):
    run_fixture(tmp_path, ["first"])
    with pytest.raises(FileExistsError):
        run_fixture(tmp_path, ["first", "new"])
    result = run_fixture(tmp_path, ["first", "new"], skip_existing=True)
    assert [job["status"] for job in result["jobs"]] == ["skipped", "success"]


def test_one_resource_receives_parallel_batch_capacity(tmp_path, monkeypatch):
    monkeypatch.setattr(os, "sched_getaffinity", lambda pid: set(range(8)), raising=False)
    result = orchestrate(
        ["first"],
        budget=Budget(10 * GIB, 6, jobs=1, worker_ram=10 * GIB),
        version="1",
        output_dir=tmp_path,
        executor=LocalExecutor(),
        poll_interval=0.02,
        worker_command=WORKER,
    )
    spec = json.loads(next(Path(result["run_dir"]).glob("first/attempt-*/spec.json")).read_text())
    assert spec["cpus"] == 6
    assert spec["allowed_cpus"] == list(range(6))
    assert spec["build"]["batch_workers"] == 6
    assert spec["build"]["resource_ram_bytes"] == 10 * GIB
    assert spec["build"]["max_batch_records"] == 50000
    assert spec["build"]["raw_batch_bytes"] == 10 * GIB // 6
    assert spec["build"]["max_batch_bytes"] == 64 * MIB


def test_external_lock_is_not_deleted(tmp_path):
    lock = tmp_path / "resources/first/.1.lock"
    lock.parent.mkdir(parents=True)
    lock.write_text("someone-else")
    with BuildLock(lock):
        result = run_fixture(tmp_path, ["first"])
    assert result["status"] == "failed"
    assert lock.read_text() == "someone-else"


def test_cancel_reaps_workers_and_leaves_no_locks(tmp_path):
    snapshots = []
    result = run_fixture(
        tmp_path,
        ["cancel", "next"],
        on_progress=snapshots.append,
        should_cancel=lambda: bool(snapshots),
    )
    assert result["status"] == "cancelled"
    assert all(j["status"] == "cancelled" for j in result["jobs"])
    assert_locks_released(tmp_path)


def test_successful_parent_cannot_leave_children_running(tmp_path):
    result = run_fixture(tmp_path, ["child"])
    pid = int(next(Path(result["run_dir"]).glob("child/attempt-*/child.pid")).read_text())
    for _ in range(20):
        try:
            if psutil.Process(pid).status() == psutil.STATUS_ZOMBIE:
                break
        except psutil.NoSuchProcess:
            break
        time.sleep(0.05)
    else:
        pytest.fail("Worker subprocess escaped its reservation")


def test_plain_console_has_progress_without_terminal_control_codes():
    stream = io.StringIO()
    console = ConsoleProgress(stream)
    console(
        {
            "elapsed_seconds": 30,
            "status": "success",
            "executor": "test",
            "budget": {"ram": GIB, "cpus": 2},
            "ram_bytes": 0,
            "reserved_bytes": 0,
            "cpu_cores": 0,
            "jobs": [{"source": "fixture", "status": "success", "message": "Done\x1b[31m"}],
        }
    )
    assert "1 done" in stream.getvalue() and "fixture" in stream.getvalue()
    assert "\x1b" not in stream.getvalue()


@pytest.mark.skipif(
    os.environ.get("OMNIPATH_TEST_SYSTEMD") != "1", reason="Explicit Linux systemd integration test"
)
def test_systemd_enforces_cpu_set_and_kernel_memory_limit(tmp_path):
    snapshots = []

    class InspectExecutor(SystemdExecutor):
        def open(self, run_id, budget):
            super().open(run_id, budget)
            props = self.ctl(
                "show",
                self.slice,
                "-p",
                "MemoryMax",
                "-p",
                "CPUQuotaPerSecUSec",
                "-p",
                "AllowedCPUs",
            ).stdout
            assert f"MemoryMax={budget.ram}" in props
            assert f"CPUQuotaPerSecUSec={budget.cpus}s" in props
            assert "AllowedCPUs=" in props

    result = orchestrate(
        ["normal", "kernel_oom"],
        budget=Budget(512 * MIB, 2, worker_ram=256 * MIB),
        version="1",
        output_dir=tmp_path,
        executor=InspectExecutor(),
        worker_command=WORKER,
        max_retries=1,
        poll_interval=0.1,
        on_progress=snapshots.append,
    )
    jobs = {j["source"]: j for j in result["jobs"]}
    assert jobs["normal"]["status"] == "success"
    assert jobs["kernel_oom"]["status"] == "failed"
    assert jobs["kernel_oom"]["attempts"] == 2
    assert jobs["kernel_oom"]["result"]["memory_failure"] is True
    assert jobs["normal"]["result"]["allowed_cpus"] == result["budget"]["allowed_cpus"]
    assert len(result["budget"]["allowed_cpus"]) < len(os.sched_getaffinity(0))
    assert_locks_released(tmp_path)


def test_cli_selects_all_resources_and_propagates_budget(monkeypatch, tmp_path):
    from omnipath_build import cli, discovery, orchestrator as module

    calls = []
    monkeypatch.setattr(
        discovery, "list_sources", lambda: [{"source": "first"}, {"source": "second"}]
    )

    def capture(sources, **kwargs):
        calls.append((sources, kwargs))
        return {"run_dir": str(tmp_path), "status": "success"}

    monkeypatch.setattr(module, "orchestrate", capture)
    assert cli.main(["run", "--all", "--version", "1", "--ram", "10GiB", "--cpus", "8"]) == 0
    sources, settings = calls[0]
    assert sources == ["first", "second"]
    assert settings["budget"].ram == 10 * GIB and settings["budget"].reserve_cpus == 1
    assert isinstance(settings["executor"], SystemdExecutor)


def test_cli_propagates_failed_run_exit_code(monkeypatch, tmp_path):
    from omnipath_build import cli, orchestrator as module

    monkeypatch.setattr(
        module, "orchestrate", lambda *a, **kw: {"run_dir": str(tmp_path), "status": "failed"}
    )
    assert cli.main(["run", "first", "--version", "1", "--ram", "1GiB", "--cpus", "2"]) == 1


def test_history_retains_working_resource_reservation(tmp_path):
    first = run_fixture(tmp_path, ["retry"])
    assert first["status"] == "success"
    learned = json.loads((tmp_path / ".build-history.json").read_text())["retry:all"]
    result = orchestrate(
        ["retry"],
        budget=Budget(GIB, 2, worker_ram=128 * MIB),
        version="2",
        output_dir=tmp_path,
        executor=LocalExecutor(),
        poll_interval=0.02,
        worker_command=WORKER,
    )
    spec = json.loads(next(Path(result["run_dir"]).glob("retry/attempt-1/spec.json")).read_text())
    assert spec["build"]["resource_ram_bytes"] == learned["reservation_bytes"]
    assert "duckdb_memory_limit" not in spec["build"]


def test_no_ineffective_memory_retry_at_full_budget_with_one_worker(tmp_path):
    result = orchestrate(
        ["retry"],
        budget=Budget(GIB, 2, worker_ram=GIB),
        version="1",
        output_dir=tmp_path,
        executor=LocalExecutor(),
        poll_interval=0.02,
        worker_command=WORKER,
    )
    (job,) = result["jobs"]
    assert job["status"] == "failed" and job["attempts"] == 1


def test_memory_retry_reduces_workers_at_full_budget(tmp_path, monkeypatch):
    monkeypatch.setattr(
        "omnipath_build.orchestrator.cpu_allocation", lambda budget: (budget, [0, 1])
    )
    result = orchestrate(
        ["retry"],
        budget=Budget(2 * GIB, 2, jobs=1, worker_ram=2 * GIB),
        version="1",
        output_dir=tmp_path,
        executor=LocalExecutor(),
        poll_interval=0.02,
        worker_command=WORKER,
    )
    (job,) = result["jobs"]
    assert job["status"] == "success" and job["attempts"] == 2
    specs = [
        json.loads(p.read_text())
        for p in sorted(Path(result["run_dir"]).glob("retry/attempt-*/spec.json"))
    ]
    assert [spec["build"]["batch_workers"] for spec in specs] == [2, 1]


def test_memory_retry_doubles_memory_and_keeps_workers(tmp_path, monkeypatch):
    monkeypatch.setattr(
        "omnipath_build.orchestrator.cpu_allocation", lambda budget: (budget, [0, 1])
    )
    result = orchestrate(
        ["retry"],
        budget=Budget(4 * GIB, 2, jobs=1, worker_ram=2 * GIB),
        version="1",
        output_dir=tmp_path,
        executor=LocalExecutor(),
        poll_interval=0.02,
        worker_command=WORKER,
    )
    (job,) = result["jobs"]
    assert job["status"] == "success" and job["attempts"] == 2
    specs = [
        json.loads(p.read_text())
        for p in sorted(Path(result["run_dir"]).glob("retry/attempt-*/spec.json"))
    ]
    assert [spec["build"]["batch_workers"] for spec in specs] == [2, 2]
    assert [spec["build"]["resource_ram_bytes"] for spec in specs] == [2 * GIB, 4 * GIB]


def test_disk_reserve_stops_workers_before_filesystem_is_full(tmp_path, monkeypatch):
    from omnipath_build import orchestrator as module

    snapshots = []
    monkeypatch.setattr(
        module.shutil, "disk_usage", lambda _: SimpleNamespace(free=0 if snapshots else GIB)
    )
    result = run_fixture(
        tmp_path, ["cancel"], on_progress=snapshots.append, min_free_disk=512 * MIB
    )
    assert result["status"] == "failed" and "Free disk" in result["error"]
    assert result["jobs"][0]["status"] == "cancelled"
    assert_locks_released(tmp_path)
