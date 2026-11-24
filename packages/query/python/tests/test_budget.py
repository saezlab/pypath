"""The defaults of the database build: threads and memory from the machine."""

from pathlib import Path

import duckdb
import pytest

from omnipath_client import ClientError, budget
from omnipath_client import database
from omnipath_client.budget import build_defaults

GB = 2**30


def machine(monkeypatch, *, cpus, memory, limit=None):
    """A machine with ``cpus`` and ``memory`` bytes available, in a container with a
    memory ``limit`` or none."""
    monkeypatch.setattr(budget, "cpu_count", lambda cgroup=None: cpus)
    monkeypatch.setattr(budget, "_system_memory", lambda: memory)
    monkeypatch.setattr(budget, "cgroup_memory_limit", lambda cgroup=None: limit)


def test_threads_are_capped_at_three(monkeypatch):
    machine(monkeypatch, cpus=64, memory=256 * GB)
    assert build_defaults()[1] == 3
    machine(monkeypatch, cpus=2, memory=256 * GB)
    assert build_defaults()[1] == 2
    machine(monkeypatch, cpus=1, memory=256 * GB)
    assert build_defaults()[1] == 1


def test_memory_is_a_share_of_the_available(monkeypatch):
    machine(monkeypatch, cpus=8, memory=16 * GB)
    assert build_defaults()[0] == f"{int(16 * 1024 * 0.6)}MiB"


def test_memory_is_never_below_the_floor(monkeypatch):
    machine(monkeypatch, cpus=4, memory=2 * GB)
    assert build_defaults()[0] == "2048MiB"


def test_a_container_limit_below_the_machine_s_memory_wins(monkeypatch):
    machine(monkeypatch, cpus=4, memory=256 * GB, limit=10 * GB)
    assert build_defaults()[0] == f"{int(10 * 1024 * 0.6)}MiB"


def test_without_a_container_limit_the_machine_counts(monkeypatch):
    machine(monkeypatch, cpus=4, memory=20 * GB)
    assert build_defaults()[0] == f"{int(20 * 1024 * 0.6)}MiB"


def test_unknown_memory_falls_back_to_a_safe_size(monkeypatch):
    machine(monkeypatch, cpus=4, memory=None)
    assert build_defaults()[0] == "4096MiB"


def test_given_settings_stay(monkeypatch):
    machine(monkeypatch, cpus=64, memory=256 * GB)
    assert build_defaults("7GB", 5) == ("7GB", 5)
    assert build_defaults("7GB") == ("7GB", 3)
    assert build_defaults(threads=1)[1] == 1
    with pytest.raises(ValueError):
        build_defaults(threads=0)


def test_cpu_quota_and_affinity(tmp_path, monkeypatch):
    monkeypatch.setattr(budget.os, "sched_getaffinity", lambda _pid: set(range(64)), raising=False)
    (tmp_path / "cpu.max").write_text("250000 100000\n")
    assert budget.cpu_count(tmp_path) == 3  # 2.5 CPUs, rounded up
    (tmp_path / "cpu.max").write_text("max 100000\n")
    assert budget.cpu_count(tmp_path) == 64
    monkeypatch.setattr(budget.os, "sched_getaffinity", lambda _pid: {0, 1}, raising=False)
    assert budget.cpu_count(tmp_path) == 2


def test_memory_limit_of_a_cgroup(tmp_path, monkeypatch):
    (tmp_path / "memory.max").write_text("max\n")
    assert budget.cgroup_memory_limit(tmp_path) is None
    (tmp_path / "memory.max").write_text(f"{8 * GB}\n")
    assert budget.cgroup_memory_limit(tmp_path) == 8 * GB
    (tmp_path / "memory.max").unlink()
    (tmp_path / "memory").mkdir()
    (tmp_path / "memory" / "memory.limit_in_bytes").write_text(f"{2**63 - 4096}\n")
    assert budget.cgroup_memory_limit(tmp_path) is None  # v1's "no limit"


def test_the_available_memory_is_the_lower_of_the_machine_and_the_container(monkeypatch):
    monkeypatch.setattr(budget, "_system_memory", lambda: 64 * GB)
    monkeypatch.setattr(budget, "cgroup_memory_limit", lambda cgroup=None: 8 * GB)
    assert budget.available_memory() == 8 * GB
    monkeypatch.setattr(budget, "cgroup_memory_limit", lambda cgroup=None: None)
    assert budget.available_memory() == 64 * GB


def test_init_builds_with_the_defaults_and_overrides(tmp_path, monkeypatch):
    machine(monkeypatch, cpus=64, memory=32 * GB)
    seen = []

    def write(target, sources, **settings):
        seen.append((settings["memory_limit"], settings["threads"]))

    monkeypatch.setattr(database, "write_database", write)
    monkeypatch.setattr(database, "_info", lambda root: None)
    installation = {"resources": [], "release": "1"}
    database.synchronize(tmp_path, installation)
    database.synchronize(tmp_path, installation, "9GB", 2)
    assert seen == [(f"{int(32 * 1024 * 0.6)}MiB", 3), ("9GB", 2)]


def test_running_out_of_memory_says_what_to_change(tmp_path, monkeypatch):
    machine(monkeypatch, cpus=64, memory=4 * GB)

    def fail(*args, **settings):
        raise duckdb.OutOfMemoryException("Out of Memory Error")

    monkeypatch.setattr(database, "write_database", fail)
    monkeypatch.setattr(database, "_info", lambda root: None)
    with pytest.raises(ClientError, match=r"memory_limit / --memory-limit.*--threads"):
        database.synchronize(tmp_path, {"resources": [], "release": "1"})
