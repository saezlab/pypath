import pytest

from omnipath_api.queries.store import budget
from omnipath_api.queries.store.budget import query_budget
from omnipath_client.budget import cgroup_memory_limit, cpu_count

SETTINGS = ("OMNIPATH_QUERY_CONCURRENCY", "OMNIPATH_DUCKDB_THREADS", "OMNIPATH_DUCKDB_MEMORY_LIMIT")


@pytest.fixture
def cgroup(tmp_path, monkeypatch):
    for name in SETTINGS:
        monkeypatch.delenv(name, raising=False)
    monkeypatch.setattr(budget.os, "sched_getaffinity", lambda pid: set(range(64)), raising=False)
    return tmp_path


def test_a_slot_of_one_thread_per_container_cpu_and_most_of_its_memory(cgroup):
    (cgroup / "cpu.max").write_text("700000 100000\n")
    (cgroup / "memory.max").write_text(f"{13 * 2**30}\n")
    b = query_budget(cgroup=cgroup)
    assert (b.slots, b.threads) == (7, 7)
    assert b.memory == f"{int(13 * 1024 * budget.MEMORY_SHARE)}MiB"


def test_a_fractional_quota_rounds_up_and_no_limit_means_the_affinity(cgroup):
    (cgroup / "cpu.max").write_text("150000 100000\n")
    assert cpu_count(cgroup) == 2
    (cgroup / "cpu.max").write_text("max 100000\n")
    assert cpu_count(cgroup) == 64
    (cgroup / "cpu.max").unlink()
    assert cpu_count(cgroup) == 64


def test_cgroup_v1_limits(cgroup):
    (cgroup / "cpu").mkdir()
    (cgroup / "cpu" / "cpu.cfs_quota_us").write_text("400000\n")
    (cgroup / "cpu" / "cpu.cfs_period_us").write_text("100000\n")
    (cgroup / "memory").mkdir()
    (cgroup / "memory" / "memory.limit_in_bytes").write_text(f"{2**63 - 4096}\n")
    assert cpu_count(cgroup) == 4
    assert cgroup_memory_limit(cgroup) is None  # v1's "no limit"
    b = query_budget(cgroup=cgroup)
    assert (b.slots, b.threads, b.memory) == (4, 4, "4GB")  # 1 GB a slot without a limit


def test_settings_override_each_part_and_empty_ones_count_as_unset(cgroup, monkeypatch):
    (cgroup / "cpu.max").write_text("700000 100000\n")
    (cgroup / "memory.max").write_text(f"{13 * 2**30}\n")
    monkeypatch.setenv("OMNIPATH_QUERY_CONCURRENCY", "3")
    monkeypatch.setenv("OMNIPATH_DUCKDB_THREADS", "2")
    monkeypatch.setenv("OMNIPATH_DUCKDB_MEMORY_LIMIT", "3500MB")
    assert query_budget(cgroup=cgroup) == budget.QueryBudget(3, 6, "10500MB")
    for name in SETTINGS:
        monkeypatch.setenv(name, " ")  # as compose passes an unset variable
    assert query_budget(cgroup=cgroup).slots == 7


def test_invalid_settings_fail(cgroup, monkeypatch):
    monkeypatch.setenv("OMNIPATH_QUERY_CONCURRENCY", "0")
    with pytest.raises(ValueError, match="concurrency"):
        query_budget(cgroup=cgroup)
    monkeypatch.setenv("OMNIPATH_QUERY_CONCURRENCY", "2")
    monkeypatch.setenv("OMNIPATH_DUCKDB_THREADS", "0")
    with pytest.raises(ValueError, match="THREADS"):
        query_budget(cgroup=cgroup)
