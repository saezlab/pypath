"""The lipid name graph is derived from the Goslin identities the release publishes."""

from contextlib import closing
import uuid

import psycopg2
from psycopg2 import sql
import pytest

from omnipath_postgres.relational import lipid_names
from omnipath_postgres.relational.db import derived_tables
from omnipath_postgres.relational.db.schema import ensure_schema


def test_a_species_identity_gives_its_class_and_totals():
    node = lipid_names.parse_identity("goslin:species:PC 34:1")
    assert (node["lipid_name"], node["lipid_level"], node["chains_listed"], node["chains_possible"]) == (
        "PC 34:1", "species", 0, 2,
    )
    assert (node["lipid_category"], node["lipid_class"], node["total_carbon"], node["total_db"]) == (
        "GP", "PC", 34, 1,
    )
    assert node["sum_formula"] == "C42H82NO8P" and node["parser_version"]


def test_a_molecular_species_lists_its_chains():
    node = lipid_names.parse_identity("goslin:molecular_species:PC 16:0_18:1")
    assert (node["lipid_level"], node["chains_listed"], node["chains_possible"]) == (
        "molecular_species", 2, 2,
    )
    assert (node["total_carbon"], node["total_db"]) == (34, 1)


def test_a_name_with_some_chains_listed_is_partially_specified():
    node = lipid_names.parse_identity("goslin:molecular_species:TG 58:12_20:0")
    assert (node["lipid_level"], node["chains_listed"], node["chains_possible"]) == (
        "partially_specified", 2, 3,
    )


def test_what_is_not_a_goslin_identity_is_not_a_node():
    for key in ("inchikey:BSYNRYMUTXBXSQ-UHFFFAOYSA-N", "goslin:species:", "goslin:nonsense", "", None):
        assert lipid_names.parse_identity(key) is None


def test_a_name_goslin_cannot_read_is_not_a_node():
    assert lipid_names.parse_identity("goslin:species:not a lipid 99:99") is None


def test_every_non_species_node_is_a_kind_of_its_species():
    nodes = {}
    for key in (
        "goslin:species:PC 34:1",
        "goslin:molecular_species:PC 16:0_18:1",
        "goslin:sn_position:PC 16:0/18:1",
        "goslin:molecular_species:PE 16:0_18:2",
    ):
        node = lipid_names.parse_identity(key)
        nodes[(node["lipid_name"], node["lipid_level"], node["chains_listed"], node["chains_possible"])] = node
    graph_nodes, edges = lipid_names.lipid_graph(nodes.values())
    species = {k for k, v in graph_nodes.items() if v["lipid_level"] == "species"}
    # PE 34:2 is no resource's name, and its species is made for the edge to attach to.
    assert ("PC 34:1", "species", 0, 2) in species and ("PE 34:2", "species", 0, 2) in species
    assert {(child[0], parent[0]) for child, parent in edges} == {
        ("PC 16:0_18:1", "PC 34:1"),
        ("PC 16:0/18:1", "PC 34:1"),
        ("PE 16:0_18:2", "PE 34:2"),
    }
    assert {e["derivation"] for e in edges.values()} == {"species_generalization"}


@pytest.fixture
def schema(postgres_dsn):
    name = "lipids_" + uuid.uuid4().hex
    with closing(psycopg2.connect(postgres_dsn)) as connection:
        try:
            ensure_schema(connection, schema=name, indexes=False, progress=False)
            with connection.cursor() as cur:
                cur.execute(
                    sql.SQL(
                        "CREATE TABLE {}.entity_reference_context (resource text, version text, "
                        "entity_key text, entity_id uuid, reference_entity_key text, "
                        "gene_reference_keys text)"
                    ).format(sql.Identifier(name))
                )
            connection.commit()
            yield connection, name
        finally:
            connection.rollback()
            with connection.cursor() as cur:
                cur.execute(sql.SQL("DROP SCHEMA IF EXISTS {} CASCADE").format(sql.Identifier(name)))
            connection.commit()


def publish(connection, name, *keys):
    with connection.cursor() as cur:
        for number, key in enumerate(keys):
            cur.execute(
                sql.SQL("INSERT INTO {}.entity_reference_context VALUES ('r', '1', %s, %s, %s, NULL)").format(
                    sql.Identifier(name)
                ),
                [f"k{number}", str(uuid.uuid4()), key],
            )
    connection.commit()


def derive(connection, name):
    with connection.cursor() as cur:
        return derived_tables._populate_lipid_identity_graph(cur, name)


def test_the_graph_is_written_from_the_published_identities(schema, monkeypatch):
    connection, name = schema

    def no_service(*args, **kwargs):
        raise AssertionError("the graph must not read a service")

    monkeypatch.setattr(derived_tables.psycopg2, "connect", no_service)
    publish(
        connection, name,
        "goslin:species:PC 34:1",
        "goslin:molecular_species:PC 16:0_18:1",
        "goslin:molecular_species:PC 16:0_18:1",  # the same identity, published by two resources
        "inchikey:BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
        "goslin:species:not a lipid 99:99",
    )
    assert derive(connection, name) == (2, 1)
    with connection.cursor() as cur:
        cur.execute(
            sql.SQL("SELECT child_name, parent_name, relation, derivation FROM {}.lipid_name_edge").format(
                sql.Identifier(name)
            )
        )
        assert cur.fetchall() == [("PC 16:0_18:1", "PC 34:1", "is_a", "species_generalization")]


def test_a_release_without_lipid_identities_gives_an_empty_graph(schema):
    connection, name = schema
    publish(connection, name, "inchikey:BSYNRYMUTXBXSQ-UHFFFAOYSA-N")
    assert derive(connection, name) == (0, 0)


def test_a_rebuild_replaces_the_graph(schema):
    connection, name = schema
    publish(connection, name, "goslin:species:PC 34:1")
    assert derive(connection, name) == (1, 0)
    assert derive(connection, name) == (1, 0)


def test_a_schema_without_the_published_references_gives_an_empty_graph(schema):
    connection, name = schema
    with connection.cursor() as cur:
        cur.execute(sql.SQL("DROP TABLE {}.entity_reference_context").format(sql.Identifier(name)))
    connection.commit()
    assert derive(connection, name) == (0, 0)
