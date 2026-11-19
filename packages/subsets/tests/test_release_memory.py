"""The memory a subset build may use comes from the environment, with a fixed default."""

import duckdb
import pytest

import metsigdb_fixture as fixture
from omnipath_subsets.release_source import open_release


def setting(memory):
    """DuckDB's own rendering of a memory limit, so the tests compare like with like."""
    con = duckdb.connect()
    con.execute("SET memory_limit = ?", [memory])
    return con.execute("SELECT current_setting('memory_limit')").fetchone()[0]


@pytest.fixture
def release(tmp_path):
    return tmp_path, fixture.write_parquet_release(tmp_path)


def limit(connection):
    return connection.execute("SELECT current_setting('memory_limit')").fetchone()[0]


def test_the_default_limit_is_four_gigabytes(release, monkeypatch):
    monkeypatch.delenv("OMNIPATH_SUBSETS_MEMORY", raising=False)
    con = open_release(*release)
    try:
        assert limit(con) == setting("4GB")
    finally:
        con.close()


def test_the_environment_sets_the_limit(release, monkeypatch):
    monkeypatch.setenv("OMNIPATH_SUBSETS_MEMORY", "24GB")
    con = open_release(*release)
    try:
        assert limit(con) == setting("24GB")
    finally:
        con.close()


def test_an_argument_wins_over_the_environment(release, monkeypatch):
    monkeypatch.setenv("OMNIPATH_SUBSETS_MEMORY", "24GB")
    con = open_release(*release, memory_limit="2GB")
    try:
        assert limit(con) == setting("2GB")
    finally:
        con.close()
