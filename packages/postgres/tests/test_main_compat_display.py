"""Namespace-local product display IDs preserve published identity and evidence."""

from importlib.resources import files
import os
import uuid

import pytest
from psycopg2 import sql
from psycopg2.extensions import parse_dsn


def test_metsigdb_chebi_display_normalization_keeps_both_original_identities():
    dsn = os.environ.get("OMNIPATH_MAIN_PARITY_DSN")
    if not dsn:
        pytest.skip("Run explicitly in the authorized private nicesrv fixture container")
    parameters = parse_dsn(dsn)
    assert parameters.get("host") == "127.0.0.1"
    assert parameters.get("port") == "5441"
    assert parameters.get("dbname") == "omnipath_migration"
    import psycopg2

    conn = psycopg2.connect(dsn)
    schema = "main_display_" + uuid.uuid4().hex
    identifier = sql.Identifier(schema)
    bare, prefixed, set_entity = [str(uuid.uuid4()) for _ in range(3)]
    created = False
    try:
        with conn.cursor() as cur:
            cur.execute("SELECT current_database()")
            assert cur.fetchone()[0] == "omnipath_migration"
            cur.execute(sql.SQL("CREATE SCHEMA {}").format(identifier))
            created = True
            # Public is excluded: unqualified publisher DDL can only target
            # this newly created fixture schema or this connection's temp tables.
            cur.execute(sql.SQL("SET search_path = {}, pg_catalog").format(identifier))
            cur.execute("""
                CREATE TABLE entity(entity_id uuid,entity_type_id int,label text,canonical_identifier text);
                CREATE TABLE vocab_entity_type(entity_type_id int,name text);
                CREATE TABLE chemical_resolution_group_member(entity_id uuid,level_id int,group_key text);
                CREATE TABLE entity_identifier_lookup(entity_id uuid,identifier_id int);
                CREATE TABLE identifier_evidence(identifier_id int,identifier_type_id int,value text);
                CREATE TABLE vocab_identifier_type(identifier_type_id int,name text);
                CREATE TABLE entity_ontology_term(term_entity_id uuid,label text);
                CREATE TABLE entity_evidence(source_id int,entity_evidence_id uuid,taxonomy_id int);
                CREATE TABLE entity_evidence_resolution(source_id int,entity_evidence_id uuid,entity_id uuid);
                CREATE TABLE metsigdb_stage(set_entity_id uuid,metabolite_entity_id uuid);
                INSERT INTO vocab_entity_type VALUES(3,'chemical_entity');
                INSERT INTO vocab_identifier_type VALUES(8,'Chebi:MI:0474');
                INSERT INTO identifier_evidence VALUES(1,8,'456'),(2,8,'CHEBI:456');
            """)
            cur.executemany(
                "INSERT INTO entity VALUES(%s,3,%s,%s)",
                [(bare, "bare label", "456"), (prefixed, "prefixed label", "CHEBI:456")],
            )
            cur.executemany(
                "INSERT INTO entity_identifier_lookup VALUES(%s,%s)", [(bare, 1), (prefixed, 2)]
            )
            cur.executemany(
                "INSERT INTO metsigdb_stage VALUES(%s,%s)",
                [(set_entity, bare), (set_entity, prefixed)],
            )
            query = (
                files("omnipath_subsets.metsigdb")
                .joinpath("sql", "publish_membership.sql")
                .read_text()
            )
            cur.execute(query, {"source_id": 9001})
            cur.execute("SELECT entity_id::text,chebi FROM metsigdb_projection")
            assert dict(cur.fetchall()) == {bare: "456", prefixed: "456"}
            cur.execute("SELECT entity_id::text,canonical_identifier FROM entity")
            assert dict(cur.fetchall()) == {bare: "456", prefixed: "CHEBI:456"}
            cur.execute(
                "SELECT identifier_id,value FROM identifier_evidence ORDER BY identifier_id"
            )
            assert cur.fetchall() == [(1, "456"), (2, "CHEBI:456")]
    finally:
        conn.rollback()
        if created:
            assert schema.startswith("main_display_") and len(schema) == 45
            with conn.cursor() as cur:
                cur.execute(sql.SQL("DROP SCHEMA IF EXISTS {} CASCADE").format(identifier))
            conn.commit()
        conn.close()
