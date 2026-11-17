"""The exclusion metabolite Parquet product over small published releases, against the rows
frozen in ``golden/`` (see its README)."""

import json
import os
from pathlib import Path

import duckdb
import pytest

from omnipath_core import subset_layout as layout
from omnipath_subsets.exclusion_metabolites import write_exclusion_metabolites
from parquet_release import write_release_from_schema
from test_exclusion_metabolites import (  # noqa: F401  (fixtures)
    cofactor_release,
    cursor_of,
    hierarchy,
    release,
    seed_directory,
)

GOLDEN = Path(__file__).parent / "golden"
MASK = "<list_version>"


def parquet_product(fixture, tmp_path, directory):
    root = tmp_path / "data"
    version = write_release_from_schema(cursor_of(fixture), root)
    layout.build_product(
        root,
        version,
        "exclusion_metabolites",
        lambda target, context: write_exclusion_metabolites(
            target, context, seed_directory=directory
        ),
    )
    product = root / "subsets" / version / "exclusion_metabolites"
    con = duckdb.connect()
    members = con.execute(
        "SELECT exclusion_class, entity_key, assignment_source, origin, match_kind, "
        "seed_namespace, seed_identifier, evidence, via_entity_key, list_version "
        f"FROM '{product}/exclusion_metabolite.parquet' ORDER BY 1, 2, 3"
    ).fetchall()
    report = con.execute(
        "SELECT kind, exclusion_class, namespace, resolved, fallback_only, unresolved, detail "
        f"FROM '{product}/exclusion_metabolite_report.parquet'"
    ).fetchall()
    return masked(members, report)


def masked(members, report):
    """Rows as plain data, with the version of the seed files (a commit or a hash) masked.

    Whether the seed files have uncommitted edits depends on the checkout that runs the
    test, so that report row is left out.
    """
    members = [list(row) for row in members]
    versions = {row[9] for row in members}
    rows = []
    for row in report:
        if row[0] == "uncommitted_inputs":
            continue
        detail = json.dumps(
            json.loads(row[6]) if isinstance(row[6], str) else row[6], sort_keys=True
        )
        for version in versions:
            detail = detail.replace(version, MASK)
        rows.append([*row[:6], detail])
    for row in members:
        row[9] = MASK
    return members, sorted(rows, key=lambda r: json.dumps(r))


def golden_of(name):
    return GOLDEN / f"exclusion_{name}.json"


@pytest.mark.parametrize("fixture_name", ["release", "hierarchy", "cofactor_release"])
@pytest.mark.parametrize("small_seeds", [False, True])
def test_the_product_equals_the_frozen_rows(fixture_name, small_seeds, request, tmp_path):
    fixture = request.getfixturevalue(fixture_name)
    directory = None
    if small_seeds:
        (tmp_path / "seeds").mkdir()
        directory = seed_directory(tmp_path / "seeds")
    members, report = parquet_product(fixture, tmp_path, directory)
    path = golden_of(f"{fixture_name}_{'small' if small_seeds else 'vendored'}")
    if os.environ.get("UPDATE_GOLDEN"):
        path.write_text(json.dumps({"members": members, "report": report}, indent=1) + "\n")
    frozen = json.loads(path.read_text())
    expected_members, expected_report = masked(frozen["members"], frozen["report"])
    assert members == expected_members
    assert report == expected_report
    if fixture_name != "release":
        assert members, "a release with a hierarchy must give rows"


def test_the_parquet_tables_carry_the_release_digest_as_their_build(release, tmp_path):  # noqa: F811
    root = tmp_path / "data"
    version = write_release_from_schema(cursor_of(release), root)
    manifest = layout.build_product(
        root, version, "exclusion_metabolites", write_exclusion_metabolites
    )
    product = root / "subsets" / version / "exclusion_metabolites"
    con = duckdb.connect()
    builds = {
        row[0]
        for row in con.execute(
            f"SELECT DISTINCT build_id FROM '{product}/exclusion_metabolite.parquet'"
        ).fetchall()
    }
    assert builds == {manifest["release_manifest_sha256"]}
