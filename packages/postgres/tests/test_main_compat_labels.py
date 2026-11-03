"""Opt-in label fixture limited to the explicitly authorized private container.

No database operation runs during collection. The root coordinator executes
this fixture only on nicesrv's separate omnipath-migration-main-postgres.
"""

import os
import uuid

import pytest
from psycopg2 import sql
from psycopg2.extensions import parse_dsn

from omnipath_postgres.relational.labels.entity_labels import populate_entity_labels


def test_published_protein_symbols_keep_main_ranking_and_identifier_fallback():
    dsn = os.environ.get("OMNIPATH_MAIN_PARITY_DSN")
    if not dsn:
        pytest.skip("Run explicitly in the private nicesrv fixture container")
    parameters = parse_dsn(dsn)
    # Guard before connecting: this cannot reach production or either pilot.
    assert parameters.get("host") == "127.0.0.1"
    assert parameters.get("port") == "5441"
    assert parameters.get("dbname") == "omnipath_migration"
    import psycopg2

    conn = psycopg2.connect(dsn)
    schema = "main_label_" + uuid.uuid4().hex
    identifier = sql.Identifier(schema)
    genes = [uuid.uuid4() for _ in range(3)]
    term, chemical = uuid.uuid4(), uuid.uuid4()
    created = False
    try:
        with conn.cursor() as cur:
            cur.execute("SELECT current_database()")
            assert cur.fetchone()[0] == "omnipath_migration"
            cur.execute(sql.SQL("CREATE SCHEMA {}").format(identifier))
            created = True
            cur.execute(
                sql.SQL("""
                CREATE TABLE {s}.vocab_entity_type (entity_type_id int, name text);
                CREATE TABLE {s}.vocab_identifier_type (identifier_type_id int, name text);
                CREATE TABLE {s}.entity (entity_id uuid,entity_type_id int,
                    canonical_identifier text,label text,label_rule text);
                CREATE TABLE {s}.identifier_evidence (identifier_id int,
                    identifier_type_id int,value text);
                CREATE TABLE {s}.entity_identifier (entity_id uuid,identifier_id int);
            """).format(s=identifier)
            )
            cur.executemany(
                sql.SQL("INSERT INTO {}.vocab_entity_type VALUES (%s,%s)").format(identifier),
                [(1, "gene"), (2, "protein"), (3, "ontology_class"), (4, "chemical_entity")],
            )
            cur.executemany(
                sql.SQL("INSERT INTO {}.vocab_identifier_type VALUES (%s,%s)").format(identifier),
                [(5, "Gene Name Primary:OM:0200"), (16, "Name:OM:0202")],
            )
            cur.executemany(
                sql.SQL("INSERT INTO {}.entity VALUES (%s,%s,%s,NULL,NULL)").format(identifier),
                [
                    (str(genes[0]), 1, "gene-id"),
                    (str(genes[1]), 2, "protein-id"),
                    (str(genes[2]), 2, "fallback-protein"),
                    (str(term), 3, "GO:0001675"),
                    (str(chemical), 4, "CHEBI:1"),
                ],
            )
            # More attestations win before shortest and alphabetical tie breaks.
            values = [
                (1, 5, "ATTESTED"),
                (2, 5, "ATTESTED"),
                (3, 5, "A"),
                (4, 5, "ZZ"),
                (5, 5, "AA"),
                (6, 5, "LONG"),
                (7, 16, "acrosome assembly"),
                (8, 16, "GO:0001675"),
                (9, 16, "a chemical name"),
            ]
            cur.executemany(
                sql.SQL("INSERT INTO {}.identifier_evidence VALUES (%s,%s,%s)").format(identifier),
                values,
            )
            cur.executemany(
                sql.SQL("INSERT INTO {}.entity_identifier VALUES (%s,%s)").format(identifier),
                [(str(genes[0]), i) for i in (1, 2, 3)]
                + [(str(genes[1]), i) for i in (4, 5, 6)]
                + [(str(term), 7), (str(term), 8), (str(chemical), 9)],
            )
        stats = populate_entity_labels(conn, schema=schema)
        with conn.cursor() as cur:
            cur.execute(
                sql.SQL(
                    "SELECT canonical_identifier,label,label_rule FROM {}.entity ORDER BY canonical_identifier"
                ).format(identifier)
            )
            # an ontology term is named; a chemical waits for its own cascade
            assert cur.fetchall() == [
                ("CHEBI:1", "CHEBI:1", "identifier_fallback"),
                ("GO:0001675", "acrosome assembly", "name"),
                ("fallback-protein", "fallback-protein", "identifier_fallback"),
                ("gene-id", "ATTESTED", "gene_symbol"),
                ("protein-id", "AA", "published_gene_symbol"),
            ]
        assert stats.gene_symbol == 1
        assert stats.published_gene_symbol == 1
        assert stats.name == 1
        assert stats.identifier_fallback == 2
        assert stats.without_label == 0
    finally:
        conn.rollback()
        if created:
            assert schema.startswith("main_label_") and len(schema) == 43
            with conn.cursor() as cur:
                cur.execute(sql.SQL("DROP SCHEMA IF EXISTS {} CASCADE").format(identifier))
            conn.commit()
        conn.close()
