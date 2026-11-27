"""The command that writes the Parquet subsets of a published release."""

import json

import pytest

from omnipath_core import subset_layout as layout
from omnipath_subsets import cli
from parquet_release import entity_row, write_release_from_schema  # noqa: F401
from omnipath_core.fixtures import write_manifest, write_resource


def publish(root, version="2026.11"):
    entities = [
        entity_row(
            "chebi:24867",
            "small_molecule",
            "chebi",
            "CHEBI:24867",
            "monoatomic ion",
            [("chebi", "CHEBI:24867", True)],
        )
    ]
    directory = root / "resources" / "chebi" / "1"
    write_resource(directory, entities, [], [])
    write_manifest(directory, "chebi", "1")
    (root / "releases").mkdir(parents=True, exist_ok=True)
    (root / "releases" / f"{version}.json").write_text(
        json.dumps({"schema_version": 1, "version": version, "resources": {"chebi": "1"}})
    )
    return version


def run(root, *args):
    return cli.main(["build", "--release", "2026.11", "--data-root", str(root), *args])


def test_the_command_writes_the_requested_product_and_its_manifest(tmp_path, capsys):
    publish(tmp_path)
    assert run(tmp_path, "--products", "exclusion_metabolites") == 0
    directory = tmp_path / "subsets" / "2026.11" / "exclusion_metabolites"
    assert (directory / "subset_manifest.json").is_file()
    assert (directory / "exclusion_metabolite.parquet").is_file()
    printed = json.loads(capsys.readouterr().out)
    assert printed["exclusion_metabolites"]["tables"]["exclusion_metabolite.parquet"]["rows"] >= 0


def test_the_command_writes_every_product_by_default(tmp_path):
    publish(tmp_path)
    assert run(tmp_path) == 0
    found = layout.present_subsets(tmp_path, "2026.11")
    assert [s.manifest["product"] for s in found] == list(layout.PRODUCTS)


def test_a_second_run_for_the_same_release_refuses_to_overwrite(tmp_path):
    publish(tmp_path)
    run(tmp_path, "--products", "exclusion_metabolites")
    before = (tmp_path / "subsets/2026.11/exclusion_metabolites/subset_manifest.json").read_bytes()
    with pytest.raises(SystemExit) as stopped:
        run(tmp_path, "--products", "exclusion_metabolites")
    assert stopped.value.code != 0
    assert (
        tmp_path / "subsets/2026.11/exclusion_metabolites/subset_manifest.json"
    ).read_bytes() == before


def test_an_unpublished_release_stops_the_command(tmp_path):
    with pytest.raises(SystemExit) as stopped:
        run(tmp_path)
    assert stopped.value.code != 0


def test_a_product_without_the_one_it_depends_on_stops_the_command(tmp_path):
    publish(tmp_path)
    with pytest.raises(SystemExit) as stopped:
        run(tmp_path, "--products", "metsigdb")
    assert stopped.value.code != 0
    assert not (tmp_path / "subsets").exists() or not any((tmp_path / "subsets").rglob("*.parquet"))


def test_the_check_command_passes_a_healthy_release_and_lists_nothing_broken(tmp_path, capsys):
    from metsigdb_fixture import write_parquet_release

    version = write_parquet_release(tmp_path)
    cli.main(["build", "--release", version, "--data-root", str(tmp_path)])
    capsys.readouterr()
    assert cli.main(["check", "--release", version, "--data-root", str(tmp_path)]) == 0
    found = json.loads(capsys.readouterr().out)
    assert set(found) == set(layout.PRODUCTS)
    assert not [f for items in found.values() for f in items if f["level"] == "error"]


def test_the_check_command_fails_for_a_release_without_subsets(tmp_path, capsys):
    publish(tmp_path)
    assert cli.main(["check", "--release", "2026.11", "--data-root", str(tmp_path)]) == 1
