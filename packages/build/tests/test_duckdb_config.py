import duckdb
import pytest

from omnipath_build.duckdb_config import build_memory_limit, configure_memory


def test_memory_budget_uses_explicit_then_environment_then_default(monkeypatch):
    monkeypatch.delenv("OMNIPATH_BUILD_DUCKDB_MEMORY", raising=False)
    assert build_memory_limit() == "1GB"
    monkeypatch.setenv("OMNIPATH_BUILD_DUCKDB_MEMORY", "256MiB")
    assert build_memory_limit() == "256MIB"
    with duckdb.connect() as con:
        configure_memory(con, "128MiB")
        assert con.execute("SELECT current_setting('memory_limit')").fetchone()[0] == "128.0 MiB"


@pytest.mark.parametrize("value", ["0GB", "-1GB", "unlimited", "1GB'; SELECT 1; --"])
def test_invalid_budget_is_rejected_before_writer_creates_files(tmp_path, value):
    from omnipath_build.writer import ParquetWriter

    with pytest.raises(ValueError, match="positive size"):
        ParquetWriter(tmp_path / "output", memory_limit=value)
    assert not (tmp_path / "output").exists()


def test_cli_forwards_memory_budget(monkeypatch, tmp_path):
    from omnipath_build import cli

    calls = []
    monkeypatch.setattr(cli, "build_resource", lambda **kwargs: calls.append(kwargs))
    assert cli.main(["build", "chembl", "--version", "1.0.0", "--memory-limit", "4GB"]) == 0
    assert calls[0]["duckdb_memory_limit"] == "4GB"


def test_writer_applies_budget_to_all_finalization_connections(tmp_path, monkeypatch):
    from omnipath_build import writer
    from test_partition_contract import run_rows, entity
    from library_fixture import build_fixture_library

    applied = []
    original = writer.configure_memory

    def record(con, limit):
        original(con, limit)
        applied.append(con.execute("SELECT current_setting('memory_limit')").fetchone()[0])

    monkeypatch.setenv("OMNIPATH_BUILD_DUCKDB_MEMORY", "128MiB")
    monkeypatch.setattr(writer, "configure_memory", record)
    library = build_fixture_library(tmp_path / "library")
    run_rows(
        tmp_path / "output",
        [("a", {"subject": entity("P04637"), "predicate": "affects", "object": entity("P02340")})],
        1,
        library,
    )
    # One disk-backed connection handles all writer stages.
    assert applied == ["128.0 MiB"]


def test_finalize_budget_leaves_headroom_under_a_hard_reservation():
    from omnipath_build.duckdb_config import finalize_budget

    gib = 1024**3
    # IntAct was killed at 4 GiB with DuckDB allowed 80% and 6 threads.
    assert finalize_budget(4 * gib, 6) == (f"{2 * gib}B", 4)
    assert finalize_budget(16 * gib, 6) == (f"{8 * gib}B", 6)
    assert finalize_budget(512 * 1024**2, 6) == (f"{256 * 1024**2}B", 1)
