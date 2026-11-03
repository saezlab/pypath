"""Shared extensions remain usable from independent main release schemas."""

import os
import uuid

import pytest
from psycopg2 import sql
from psycopg2.extensions import parse_dsn

from omnipath_postgres.relational.db.bitmaps import _create_bitmap_tables
from omnipath_postgres.relational.db.extensions import ensure_public_extension


class ExtensionCursor:
    def __init__(self, namespace):
        self.namespace = namespace
        self.queries = []

    def execute(self, query, parameters=None):
        self.queries.append((query, parameters))

    def fetchone(self):
        return (self.namespace,)


def test_extension_creation_targets_public_and_preserves_existing_install():
    cursor = ExtensionCursor("public")
    ensure_public_extension(cursor, "roaringbitmap")
    assert str(cursor.queries[0][0]) == str(
        sql.SQL("CREATE EXTENSION IF NOT EXISTS {} WITH SCHEMA public").format(
            sql.Identifier("roaringbitmap")
        )
    )
    assert cursor.queries[1][1] == ("roaringbitmap",)
    assert not any(
        "ALTER EXTENSION" in str(query) or "DROP EXTENSION" in str(query)
        for query, _params in cursor.queries
    )


def test_existing_wrong_namespace_requires_explicit_setup():
    cursor = ExtensionCursor("earlier_release")
    with pytest.raises(ValueError, match="Relocate it explicitly"):
        ensure_public_extension(cursor, "roaringbitmap")
    assert len(cursor.queries) == 2


def test_two_release_schemas_reuse_public_roaringbitmap_type():
    dsn = os.environ.get("OMNIPATH_MAIN_PARITY_DSN")
    if not dsn:
        pytest.skip("Run explicitly in the authorized private nicesrv fixture container")
    parameters = parse_dsn(dsn)
    assert parameters.get("host") == "127.0.0.1"
    assert parameters.get("port") == "5441"
    assert parameters.get("dbname") == "omnipath_migration"
    import psycopg2

    conn = psycopg2.connect(dsn)
    schemas = ["main_bitmap_" + uuid.uuid4().hex for _ in range(2)]
    try:
        with conn.cursor() as cursor:
            cursor.execute("SELECT current_database()")
            assert cursor.fetchone()[0] == "omnipath_migration"
            cursor.execute(
                "SELECT n.nspname FROM pg_extension e JOIN pg_namespace n "
                "ON n.oid=e.extnamespace WHERE e.extname='roaringbitmap'"
            )
            # This fixture never installs/moves an extension or uses old schemas.
            assert cursor.fetchone() == ("public",)
            for schema in schemas:
                identifier = sql.Identifier(schema)
                cursor.execute(sql.SQL("CREATE SCHEMA {}").format(identifier))
                cursor.execute(sql.SQL("SET search_path = {}, public").format(identifier))
                _create_bitmap_tables(cursor, schema)
                cursor.execute(
                    sql.SQL("""INSERT INTO {}.facet_entity_bitmap
                    VALUES ('fixture','value',rb_build(ARRAY[1,2,3]),3)""").format(identifier)
                )
                cursor.execute(
                    sql.SQL(
                        "SELECT rb_cardinality(entity_bitmap) FROM {}.facet_entity_bitmap"
                    ).format(identifier)
                )
                assert cursor.fetchone() == (3,)
                cursor.execute(
                    """SELECT tn.nspname, t.typname
                    FROM pg_attribute a JOIN pg_class c ON c.oid=a.attrelid
                    JOIN pg_namespace cn ON cn.oid=c.relnamespace
                    JOIN pg_type t ON t.oid=a.atttypid
                    JOIN pg_namespace tn ON tn.oid=t.typnamespace
                    WHERE cn.nspname=%s AND c.relname='facet_entity_bitmap'
                      AND a.attname='entity_bitmap' """,
                    (schema,),
                )
                assert cursor.fetchone() == ("public", "roaringbitmap")
    finally:
        # All fixture DDL/data is uncommitted and rolls back together.
        conn.rollback()
        conn.close()
