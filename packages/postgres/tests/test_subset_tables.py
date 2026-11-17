"""The tables the subsets are loaded into: their shape, vocabularies and review view."""

import uuid

import psycopg2
from psycopg2 import sql
import pytest

from omnipath_postgres.subset_loading import create_tables

CLASSES = ("ions", "atoms", "small_molecules", "xenobiotics", "cofactors", "cofactors_extended")
WATER, PROTON = str(uuid.UUID(int=1)), str(uuid.UUID(int=2))

ROW = (
    "INSERT INTO exclusion_metabolite(entity_id, exclusion_class, assignment_source, origin, "
    "match_kind, evidence, via_entity_id, list_version, build_id) "
    "VALUES (%s, %s, %s, %s, %s, 'e', %s, 'v', 'b')"
)


def cursor_of(release):
    connection, schema, _ = release
    create_tables(connection, schema, "exclusion_metabolites")
    cur = connection.cursor()
    cur.execute(sql.SQL("SET search_path={},public").format(sql.Identifier(schema)))
    return cur


def columns(cur, table):
    cur.execute(
        "SELECT column_name, data_type, is_nullable FROM information_schema.columns "
        "WHERE table_schema = current_schema() AND table_name = %s ORDER BY ordinal_position",
        [table],
    )
    return cur.fetchall()


def test_the_product_creates_both_tables_with_the_contracted_shape(release):
    cur = cursor_of(release)
    assert columns(cur, "exclusion_metabolite") == [
        ("entity_id", "uuid", "NO"),
        ("exclusion_class", "text", "NO"),
        ("assignment_source", "text", "NO"),
        ("origin", "text", "NO"),
        ("match_kind", "text", "NO"),
        ("seed_namespace", "text", "YES"),
        ("seed_identifier", "text", "YES"),
        ("evidence", "text", "NO"),
        ("via_entity_id", "uuid", "YES"),
        ("list_version", "text", "NO"),
        ("build_id", "text", "NO"),
    ]
    assert columns(cur, "exclusion_metabolite_report") == [
        ("kind", "text", "NO"),
        ("exclusion_class", "text", "YES"),
        ("namespace", "text", "YES"),
        ("resolved", "integer", "YES"),
        ("fallback_only", "integer", "YES"),
        ("unresolved", "integer", "YES"),
        ("detail", "jsonb", "YES"),
        ("build_id", "text", "NO"),
    ]
    cur.execute(
        "SELECT a.attname FROM pg_index i JOIN pg_attribute a ON a.attrelid = i.indrelid "
        "AND a.attnum = ANY(i.indkey) WHERE i.indrelid = 'exclusion_metabolite'::regclass "
        "AND i.indisprimary ORDER BY array_position(i.indkey::int[], a.attnum::int)"
    )
    assert [row[0] for row in cur.fetchall()] == [
        "entity_id",
        "exclusion_class",
        "assignment_source",
    ]
    cur.execute(
        "SELECT pg_get_indexdef(indexrelid) FROM pg_index "
        "WHERE indrelid = 'exclusion_metabolite'::regclass AND NOT indisprimary"
    )
    definitions = [row[0] for row in cur.fetchall()]
    assert any("(exclusion_class, entity_id)" in definition for definition in definitions)


@pytest.mark.parametrize(
    "row",
    [
        (WATER, "lipids", "metaproviz_seed", "curated", "entity", None),
        (WATER, "ions", "somewhere", "curated", "entity", None),
        (WATER, "ions", "metaproviz_seed", "invented", "entity", None),
        (WATER, "ions", "metaproviz_seed", "curated", "fuzzy", None),
        (WATER, "ions", "metaproviz_seed", "curated", "structure_key", None),
    ],
)
def test_the_table_refuses_values_outside_its_vocabularies(release, row):
    cur = cursor_of(release)
    with pytest.raises(psycopg2.errors.CheckViolation):
        cur.execute(ROW, row)


def test_the_table_accepts_every_class_and_a_key_row_with_its_member(release):
    cur = cursor_of(release)
    for exclusion_class in CLASSES:
        cur.execute(ROW, (WATER, exclusion_class, "metaproviz_seed", "curated", "entity", None))
    cur.execute(ROW, (PROTON, "ions", "chebi_hierarchy", "derived", "structure_key", WATER))
    with pytest.raises(psycopg2.errors.UniqueViolation):
        cur.execute(ROW, (WATER, "ions", "metaproviz_seed", "curated", "entity", None))


def test_applying_the_table_definition_twice_changes_nothing(release):
    cur = cursor_of(release)
    cur.execute(ROW, (WATER, "ions", "metaproviz_seed", "curated", "entity", None))
    release[0].commit()
    before = columns(cur, "exclusion_metabolite")
    cur = cursor_of(release)
    assert columns(cur, "exclusion_metabolite") == before
    cur.execute("SELECT count(*) FROM exclusion_metabolite")
    assert cur.fetchone()[0] == 1


def test_the_review_view_lists_members_with_their_label_and_identifiers(release):
    cur = cursor_of(release)
    cur.execute("INSERT INTO entity VALUES (%s, 'water', 'CHEBI:15377')", [WATER])
    cur.execute(
        "INSERT INTO vocab_identifier_type VALUES (1, 'Chebi:MI:0474'), (2, 'Kegg Compound:MI:2012')"
    )
    for number, (type_id, value) in enumerate([(1, "CHEBI:15377"), (2, "C00001")], start=1):
        identifier = str(uuid.UUID(int=1000 + number))
        cur.execute(
            "INSERT INTO identifier_evidence VALUES (%s, %s, %s)", [identifier, type_id, value]
        )
        cur.execute("INSERT INTO entity_identifier_lookup VALUES (%s, %s)", [WATER, identifier])
    cur.execute(ROW, (WATER, "small_molecules", "metaproviz_seed", "curated", "entity", None))
    cur.execute(
        "SELECT exclusion_class, label, chebi, kegg, assignment_source, origin, match_kind "
        "FROM exclusion_metabolite_members"
    )
    assert cur.fetchall() == [
        ("small_molecules", "water", "15377", "C00001", "metaproviz_seed", "curated", "entity")
    ]


def test_the_membership_table_refuses_a_set_type_or_resource_outside_the_contract(release):
    connection, schema, _ = release
    create_tables(connection, schema, "metsigdb")
    cur = connection.cursor()
    cur.execute(sql.SQL("SET search_path={},public").format(sql.Identifier(schema)))
    insert = (
        "INSERT INTO metsigdb_membership(resource, set_source_id, metabolite_entity_id, "
        "metabolite_label, metabolite_entity_type, has_structure, exclusion_classes, set_type, "
        "set_size, provenance_source, build_id) "
        "VALUES (%s, 's', %s, 'm', 'small_molecule', false, '{}', %s, 1, 'p', 'b')"
    )
    cur.execute(insert, ("Reactome", WATER, "pathway"))
    for row in (("Reactome", WATER, "tissue"), ("SMPDB", PROTON, "pathway")):
        with pytest.raises(psycopg2.errors.CheckViolation):
            cur.execute(insert, row)
        connection.rollback()
        cur.execute(sql.SQL("SET search_path={},public").format(sql.Identifier(schema)))
