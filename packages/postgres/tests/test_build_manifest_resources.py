"""The build manifest names each resource's version and input-module commit from the release."""

from contextlib import closing
import json
import uuid

import psycopg2
from psycopg2 import sql
import pytest

from omnipath_postgres.relational.db import resources


@pytest.fixture
def schema(postgres_dsn):
    name = "manifest_" + uuid.uuid4().hex
    with closing(psycopg2.connect(postgres_dsn)) as connection:
        try:
            with connection.cursor() as cur:
                cur.execute(sql.SQL("CREATE SCHEMA {}").format(sql.Identifier(name)))
                cur.execute(sql.SQL("SET search_path={},public").format(sql.Identifier(name)))
                cur.execute(
                    """
                    CREATE TABLE resources(resource_id text PRIMARY KEY, entity_count int,
                      interaction_count int, association_count int, identifier_count int,
                      ontology_term_count int, input_module_commit text, input_module_dirty boolean);
                    CREATE TABLE parquet_resource(resource text PRIMARY KEY, version text,
                      manifest_sha256 text, manifest jsonb);
                    CREATE TABLE build_capability(capability text PRIMARY KEY, available boolean,
                      provider text, reason text);
                    """
                )
            connection.commit()
            yield connection, name
        finally:
            connection.rollback()
            with connection.cursor() as cur:
                cur.execute(sql.SQL("DROP SCHEMA {} CASCADE").format(sql.Identifier(name)))
            connection.commit()


def add(connection, schema, name, *, checkout=("checkout", False), release=None):
    with connection.cursor() as cur:
        cur.execute(sql.SQL("SET search_path={},public").format(sql.Identifier(schema)))
        cur.execute(
            "INSERT INTO resources VALUES (%s, 1, 2, 3, 4, 5, %s, %s)", [name, *checkout]
        )
        if release is not None:
            version, manifest = release
            cur.execute(
                "INSERT INTO parquet_resource VALUES (%s, %s, 'digest', %s)",
                [name, version, json.dumps(manifest)],
            )
    connection.commit()


def manifest_resources(connection, schema):
    with connection.cursor() as cur:
        cur.execute(sql.SQL("SET search_path={},public").format(sql.Identifier(schema)))
        return resources._manifest_resources(cur, schema)


def software(commit, dirty=False):
    return {
        "provenance": {
            "software": {"packages": {"omnipath-resources": {"git": {"commit": commit, "dirty": dirty}}}}
        }
    }


def test_a_resource_of_the_release_has_its_pinned_version_and_the_commit_it_was_built_at(schema):
    connection, name = schema
    add(connection, name, "signor", release=("2.1.0", software("release-commit", True)))
    (found,) = manifest_resources(connection, name)
    assert found == {
        "name": "signor",
        "record_count": 15,
        "version": "2.1.0",
        "input_module_commit": "release-commit",
        "input_module_dirty": True,
    }


def test_a_resource_the_release_does_not_pin_keeps_the_checkouts_commit_and_no_version(schema):
    connection, name = schema
    add(connection, name, "orphan", checkout=("checkout-commit", True))
    (found,) = manifest_resources(connection, name)
    assert (found["version"], found["input_module_commit"], found["input_module_dirty"]) == (
        None,
        "checkout-commit",
        True,
    )


def test_a_build_that_recorded_no_revision_keeps_the_checkouts_commit(schema):
    connection, name = schema
    add(connection, name, "signor", checkout=("checkout-commit", False), release=("1", {}))
    (found,) = manifest_resources(connection, name)
    assert (found["version"], found["input_module_commit"]) == ("1", "checkout-commit")


def test_resources_come_in_name_order(schema):
    connection, name = schema
    for resource in ("zeta", "alpha", "mid"):
        add(connection, name, resource, release=("1", software(resource)))
    assert [r["name"] for r in manifest_resources(connection, name)] == ["alpha", "mid", "zeta"]


def test_a_schema_without_the_release_tables_still_lists_its_resources(schema):
    connection, name = schema
    add(connection, name, "signor")
    with connection.cursor() as cur:
        cur.execute(sql.SQL("DROP TABLE {}.parquet_resource").format(sql.Identifier(name)))
    connection.commit()
    (found,) = manifest_resources(connection, name)
    assert found["version"] is None and found["input_module_commit"] == "checkout"


def emitted(connection, schema):
    with connection.cursor() as cur:
        cur.execute(sql.SQL("SET search_path={},public").format(sql.Identifier(schema)))
        resources.emit_build_manifest(connection, schema=schema)
        cur.execute(
            sql.SQL("SELECT * FROM {}.build_manifest").format(sql.Identifier(schema))
        )
        names = [column.name for column in cur.description]
        return dict(zip(names, cur.fetchone()))


def test_the_manifest_names_no_translation_tables_the_build_no_longer_uses(schema):
    connection, name = schema
    add(connection, name, "signor", release=("2.1.0", software("release-commit")))
    row = emitted(connection, name)
    assert "translation_tables" not in row
    assert row["resources"][0]["version"] == "2.1.0"
    assert row["resources"][0]["input_module_commit"] == "release-commit"


def test_the_build_id_follows_the_pinned_release_not_the_checkout(schema):
    connection, name = schema
    add(connection, name, "signor", checkout=("one", False), release=("1", software("pinned")))
    first = emitted(connection, name)["build_id"]
    with connection.cursor() as cur:
        cur.execute(sql.SQL("SET search_path={},public").format(sql.Identifier(name)))
        cur.execute("UPDATE resources SET input_module_commit = 'another checkout'")
    connection.commit()
    assert emitted(connection, name)["build_id"] == first


def built_with(**dependencies):
    return {"provenance": {"dependencies": dependencies}}


def capabilities(connection, schema):
    with connection.cursor() as cur:
        cur.execute(sql.SQL("SET search_path={},public").format(sql.Identifier(schema)))
        return {c["capability"]: c for c in resources._build_capabilities(cur, schema)}


def test_lipid_naming_and_structure_keys_follow_what_the_pinned_builds_recorded(schema):
    connection, name = schema
    add(connection, name, "swisslipids", release=("1", built_with(pygoslin="2.0.1", rdkit="2025.3")))
    found = capabilities(connection, name)
    for capability in ("lipid_nomenclature", "structure_key_computation"):
        assert found[capability]["available"] is True
        assert found[capability]["reason"] is None
        assert "release" in found[capability]["provider"]


def test_a_toolkit_no_pinned_build_recorded_is_unavailable_with_the_reason(schema):
    connection, name = schema
    add(connection, name, "signor", release=("1", built_with(pygoslin=None, rdkit="2025.3")))
    found = capabilities(connection, name)
    assert found["structure_key_computation"]["available"] is True
    assert found["lipid_nomenclature"]["available"] is False
    assert "pygoslin" in found["lipid_nomenclature"]["reason"]


def test_a_schema_without_the_release_says_what_is_missing(schema):
    connection, name = schema
    with connection.cursor() as cur:
        cur.execute(sql.SQL("DROP TABLE {}.parquet_resource").format(sql.Identifier(name)))
    connection.commit()
    found = capabilities(connection, name)
    assert found["lipid_nomenclature"]["available"] is False
    assert "parquet_resource" in found["lipid_nomenclature"]["reason"]


def test_the_capabilities_do_not_read_the_utils_database():
    import inspect

    assert "utils_db_url" not in inspect.signature(resources._build_capabilities).parameters
    assert "utils_db_url" not in inspect.signature(resources.emit_build_manifest).parameters
    assert "build_info" not in inspect.getsource(resources)


def test_the_manifest_stores_the_capabilities_it_found(schema):
    connection, name = schema
    add(connection, name, "swisslipids", release=("1", built_with(pygoslin="2.0.1", rdkit="2025.3")))
    emitted(connection, name)
    with connection.cursor() as cur:
        cur.execute(
            sql.SQL("SELECT capability, available FROM {}.build_capability ORDER BY 1").format(
                sql.Identifier(name)
            )
        )
        stored = dict(cur.fetchall())
    assert stored["lipid_nomenclature"] and stored["structure_key_computation"]
