"""Checkpoint and finish options cannot silently ignore load/audit arguments."""

import pytest

from omnipath_postgres import cli


@pytest.mark.parametrize(
    "extra",
    [
        ["manifest.json"],
        ["--validate-source-records"],
        ["--batch-size", "1"],
        ["--checkpoint-base"],
        ["--defer-constraints"],
    ],
)
def test_finish_rejects_inapplicable_options_without_connecting(monkeypatch, extra):
    def unexpected(*args, **kwargs):
        pytest.fail("Invalid CLI selection must not connect")

    monkeypatch.setattr(cli, "finish_release", unexpected)
    monkeypatch.setattr(cli, "load_release", unexpected)
    with pytest.raises(SystemExit) as error:
        cli.main(["--finish", "--database-url", "unused", *extra])
    assert error.value.code == 2


def test_cli_forwards_bulk_load_and_checkpoint_options(monkeypatch, capsys):
    from omnipath_postgres import cli
    from omnipath_postgres.loader import LoadResult

    calls = []

    def load(*args, **kwargs):
        calls.append(kwargs)
        return LoadResult("test", "release", "hash", {}, {}, {}, (), {})

    monkeypatch.setattr(cli, "load_release", load)
    assert (
        cli.main(
            [
                "release.json",
                "--data-root",
                "data",
                "--database-url",
                "unused",
                "--defer-constraints",
                "--checkpoint-base",
            ]
        )
        == 0
    )
    assert calls[0]["defer_constraints"] is True
    assert calls[0]["base_only"] is False
    assert tuple(calls[0]["products"]) == (
        "exclusion_metabolites",
        "metsigdb",
        "network_views",
        "cosmos",
    )
    capsys.readouterr()


def test_finish_accepts_the_exclusion_metabolites_product(monkeypatch, capsys):
    from omnipath_postgres.loader import LoadResult

    calls = []

    def finish(database_url, **kwargs):
        calls.append(kwargs)
        return LoadResult("test", "release", "hash", {}, {}, {}, (), {})

    monkeypatch.setattr(cli, "finish_release", finish)
    assert (
        cli.main(["--finish", "--database-url", "unused", "--products", "exclusion_metabolites"])
        == 0
    )
    assert tuple(calls[0]["products"]) == ("exclusion_metabolites",)
    capsys.readouterr()


def test_finish_forwards_the_data_root_that_holds_the_subsets(monkeypatch, capsys):
    from omnipath_postgres.loader import LoadResult

    calls = []

    def finish(database_url, **kwargs):
        calls.append(kwargs)
        return LoadResult("test", "release", "hash", {}, {}, {}, (), {})

    monkeypatch.setattr(cli, "finish_release", finish)
    assert (
        cli.main(
            [
                "--finish",
                "--database-url",
                "unused",
                "--data-root",
                "data",
                "--duckdb-threads",
                "2",
                "--memory-limit",
                "2GB",
                "--temp-directory",
                "tmp",
            ]
        )
        == 0
    )
    assert calls[0]["subsets_root"] == "data"
    assert (calls[0]["duckdb_threads"], calls[0]["memory_limit"]) == (2, "2GB")
    assert calls[0]["temp_directory"] == "tmp"
    capsys.readouterr()
