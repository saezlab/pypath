"""Release checks over the subsets of a release: what must hold of any build, relative to it.

A healthy release passes; each check has a case that breaks exactly it. No count is frozen:
the checks compare the tables with each other and with the release they came from.
"""

import json
from pathlib import Path
import shutil

import duckdb
import pytest

from omnipath_subsets import checks
from omnipath_core.subset_layout import PRODUCTS
from omnipath_subsets.parquet_products import build_release_subsets
from omnipath_subsets.release_source import open_release
import metsigdb_fixture as fixture

MEMBERSHIP = "metsigdb_membership.parquet"
MEMBERS = "exclusion_metabolite.parquet"
REPORT = "exclusion_metabolite_report.parquet"


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    root = tmp_path_factory.mktemp("release")
    version = fixture.write_parquet_release(root)
    build_release_subsets(root, version)
    return root, version


@pytest.fixture
def release_connection(built):
    root, version = built
    con = open_release(root, version)
    yield con
    con.close()


def digest(built):
    root, version = built
    return json.loads(
        (root / "subsets" / version / "metsigdb" / "subset_manifest.json").read_text()
    )["release_manifest_sha256"]


def corrupted(built, tmp_path, table, statement, product="metsigdb"):
    """A copy of one subset table with ``statement`` (SQL over ``t``) applied."""
    root, version = built
    source = root / "subsets" / version / product / table
    target = tmp_path / table
    shutil.copy(source, target)
    con = duckdb.connect()
    con.execute(f"CREATE TABLE t AS SELECT * FROM read_parquet('{target}')")
    for part in statement.split(";"):
        if part.strip():
            con.execute(part)
    con.execute(f"COPY t TO '{target}' (FORMAT parquet)")
    return target


def errors(findings):
    return {f.check for f in findings if f.level == "error"}


# ----------------------------------------------------------------- a healthy release


def test_a_healthy_release_has_no_error(built):
    root, version = built
    found = checks.check_release(root, version)
    assert set(found) == set(PRODUCTS)
    assert {name: errors(items) for name, items in found.items()} == {
        name: set() for name in PRODUCTS
    }


def test_unnamed_sets_of_the_resources_that_may_lack_names_are_reported_not_refused(built):
    root, version = built
    kegg = [f for f in checks.check_release(root, version)["metsigdb"] if f.check == "set_labels"]
    assert kegg and all(f.level == "warning" for f in kegg)


# ----------------------------------------------------------------- MetSigDB


def run(built, table, release_connection=None, **kwargs):
    return checks.check_metsigdb(
        table, release=release_connection, build_id=digest(built), **kwargs
    )


def test_a_membership_in_the_table_twice_is_refused(built, tmp_path):
    table = corrupted(built, tmp_path, MEMBERSHIP, "INSERT INTO t SELECT * FROM t LIMIT 1")
    assert "no_duplicates" in errors(run(built, table))


@pytest.mark.parametrize("column", ["metabolite_label", "set_type", "provenance_source", "build_id"])
def test_a_row_without_a_mandatory_field_is_refused(built, tmp_path, column):
    table = corrupted(built, tmp_path, MEMBERSHIP, f"UPDATE t SET {column} = NULL WHERE resource = 'Reactome'")
    assert "mandatory_fields" in errors(run(built, table))


def test_a_set_type_that_is_not_the_resources_is_refused(built, tmp_path):
    table = corrupted(built, tmp_path, MEMBERSHIP, "UPDATE t SET set_type = 'disease' WHERE resource = 'Reactome'")
    assert "set_types" in errors(run(built, table))


def test_a_set_size_that_is_not_the_sets_population_is_refused(built, tmp_path):
    table = corrupted(built, tmp_path, MEMBERSHIP, "UPDATE t SET set_size = set_size + 1 WHERE resource = 'Reactome'")
    assert "set_size" in errors(run(built, table))


def test_a_build_stamp_of_another_release_is_refused(built, tmp_path):
    table = corrupted(built, tmp_path, MEMBERSHIP, "UPDATE t SET build_id = 'x'")
    assert "build_stamp" in errors(run(built, table))


def test_a_set_that_reactome_does_not_name_is_refused(built, tmp_path):
    table = corrupted(built, tmp_path, MEMBERSHIP, "UPDATE t SET set_label = NULL WHERE resource = 'Reactome'")
    assert "set_labels" in errors(run(built, table))


def test_a_classyfire_set_that_is_not_a_chemont_term_is_refused(built, tmp_path):
    table = corrupted(built, tmp_path, MEMBERSHIP, "UPDATE t SET set_source_id = 'OMIM:1' WHERE resource = 'ClassyFire' AND set_source_id = 'CHEMONTID:L'")
    assert "chemont_prefix" in errors(run(built, table))


def test_a_macdb_set_id_without_its_prefix_is_refused(built, tmp_path):
    table = corrupted(built, tmp_path, MEMBERSHIP, "UPDATE t SET set_source_id = replace(set_source_id, 'MACDB:', '') WHERE resource = 'MACdb'")
    assert "macdb_prefix" in errors(run(built, table))


def test_an_organism_on_a_kegg_reference_map_is_refused(built, tmp_path):
    table = corrupted(built, tmp_path, MEMBERSHIP, "UPDATE t SET organism = 10090 WHERE set_source_id = 'rn00010'")
    assert "organism" in errors(run(built, table))


def test_an_overview_map_that_is_not_marked_as_one_is_refused(built, tmp_path):
    table = corrupted(built, tmp_path, MEMBERSHIP, "UPDATE t SET set_sub_type = 'metabolic_map' WHERE set_source_id = 'rn01100'")
    assert "overview_maps" in errors(run(built, table))


def test_a_pathway_marked_as_an_overview_map_is_refused(built, tmp_path):
    table = corrupted(built, tmp_path, MEMBERSHIP, "UPDATE t SET set_sub_type = 'overview_map' WHERE set_source_id = 'rn00010'")
    assert "overview_maps" in errors(run(built, table))


def test_an_unknown_exclusion_class_is_refused(built, tmp_path):
    table = corrupted(built, tmp_path, MEMBERSHIP, "UPDATE t SET exclusion_classes = ['lipids']")
    assert "exclusion_classes" in errors(run(built, table))


def test_a_structure_flag_that_contradicts_the_identifiers_is_refused(built, tmp_path):
    table = corrupted(built, tmp_path, MEMBERSHIP, "UPDATE t SET has_structure = NOT has_structure")
    assert "has_structure" in errors(run(built, table))


def test_a_metabolite_the_release_does_not_have_is_refused(built, tmp_path, release_connection):
    table = corrupted(built, tmp_path, MEMBERSHIP, "UPDATE t SET metabolite_entity_key = 'nowhere' WHERE metabolite_entity_key = 'chemical_b'")
    assert "canonical_entity" in errors(run(built, table, release_connection))


def test_a_member_that_is_not_a_chemical_is_refused(built, tmp_path, release_connection):
    table = corrupted(built, tmp_path, MEMBERSHIP, "UPDATE t SET metabolite_entity_key = 'protein' WHERE metabolite_entity_key = 'chemical_b'")
    assert "canonical_entity" in errors(run(built, table, release_connection))


def test_the_entity_check_is_left_out_without_the_release(built, tmp_path):
    table = corrupted(built, tmp_path, MEMBERSHIP, "UPDATE t SET metabolite_entity_key = 'nowhere'")
    assert "canonical_entity" not in {f.check for f in run(built, table)}


# ----------------------------------------------------------------- exclusion metabolites


def run_exclusion(built, members, report, release_connection=None):
    return checks.check_exclusion(
        members, report, release=release_connection, build_id=digest(built)
    )


@pytest.fixture
def exclusion(built, tmp_path):
    """The exclusion tables of the release with a member in every class, as copies."""
    members, report = tmp_path / MEMBERS, tmp_path / REPORT
    con = duckdb.connect()
    con.execute(
        "CREATE TABLE m AS SELECT 'chemical_a' AS entity_key, c AS exclusion_class, "
        "'project_seed' AS assignment_source, 'curated' AS origin, 'entity' AS match_kind, "
        "'chebi' AS seed_namespace, 'CHEBI:1' AS seed_identifier, 'seed' AS evidence, "
        "NULL::VARCHAR AS via_entity_key, 'v' AS list_version, ? AS build_id "
        "FROM unnest(['ions','atoms','small_molecules','xenobiotics','cofactors','cofactors_extended']) AS u(c)",
        [digest(built)],
    )
    con.execute(f"COPY m TO '{members}' (FORMAT parquet)")
    # The report agrees with the members: one entity in each class.
    con.execute(
        "CREATE TABLE r AS SELECT 'class_size' AS kind, exclusion_class, NULL::VARCHAR AS namespace, "
        "1 AS resolved, NULL::INTEGER AS fallback_only, NULL::INTEGER AS unresolved, "
        "'{}' AS detail, build_id FROM m"
    )
    con.execute(f"COPY r TO '{report}' (FORMAT parquet)")
    return members, report, con


def alter(con, path, table, statement):
    for part in statement.split(";"):
        if part.strip():
            con.execute(part)
    con.execute(f"COPY {table} TO '{path}' (FORMAT parquet)")


def test_healthy_exclusion_tables_pass(built, exclusion):
    members, report, _ = exclusion
    assert errors(run_exclusion(built, members, report)) == set()


def test_a_class_outside_the_vocabulary_is_refused(built, exclusion):
    members, report, con = exclusion
    alter(con, members, "m", "UPDATE m SET exclusion_class = 'lipids' WHERE exclusion_class = 'ions'")
    assert "classes" in errors(run_exclusion(built, members, report))


def test_a_member_without_a_class_is_refused(built, exclusion):
    members, report, con = exclusion
    alter(con, members, "m", "UPDATE m SET exclusion_class = NULL WHERE exclusion_class = 'ions'")
    assert "classes" in errors(run_exclusion(built, members, report))


def test_an_entity_listed_twice_in_a_class_by_one_source_is_refused(built, exclusion):
    members, report, con = exclusion
    alter(con, members, "m", "INSERT INTO m SELECT * FROM m LIMIT 1")
    assert "no_duplicates" in errors(run_exclusion(built, members, report))


def test_a_structure_key_row_without_its_member_is_refused(built, exclusion):
    members, report, con = exclusion
    alter(con, members, "m", "UPDATE m SET match_kind = 'structure_key' WHERE exclusion_class = 'ions'")
    assert "structure_key_members" in errors(run_exclusion(built, members, report))


def test_a_class_without_a_size_row_in_the_report_is_refused(built, exclusion):
    members, report, con = exclusion
    alter(con, report, "r", "DELETE FROM r WHERE kind = 'class_size' AND exclusion_class = 'ions'")
    assert "report_classes" in errors(run_exclusion(built, members, report))


def test_a_class_size_the_members_do_not_have_is_refused(built, exclusion):
    members, report, con = exclusion
    alter(con, report, "r", "UPDATE r SET resolved = resolved + 5 WHERE kind = 'class_size' AND exclusion_class = 'ions'")
    assert "report_classes" in errors(run_exclusion(built, members, report))


def test_a_build_stamp_of_another_release_on_the_exclusion_tables_is_refused(built, exclusion):
    members, report, con = exclusion
    alter(con, members, "m", "UPDATE m SET build_id = 'x'")
    assert "build_stamp" in errors(run_exclusion(built, members, report))


def test_a_member_the_release_does_not_have_is_refused(built, exclusion, release_connection):
    members, report, con = exclusion
    alter(con, members, "m", "UPDATE m SET entity_key = 'nowhere' WHERE exclusion_class = 'ions'")
    assert "canonical_entity" in errors(run_exclusion(built, members, report, release_connection))
