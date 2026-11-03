"""Partitioned late indexes must retain main's exact valid parent/child contract."""

from collections import Counter
from contextlib import closing
import uuid

import psycopg2
from psycopg2 import sql
import pytest

from omnipath_postgres import loader as aligned_loader


INDEX = "entity_evidence_resolution_molecular_type_idx"
TABLE = "entity_evidence_resolution"


@pytest.fixture
def index_database(postgres_dsn):
    """Four declared rows, two source partitions and a default partition."""
    schema = "index_replay_" + uuid.uuid4().hex
    with closing(psycopg2.connect(postgres_dsn)) as connection:

        def query(text):
            return sql.SQL(text).format(s=sql.Identifier(schema))

        with connection.cursor() as cursor:
            cursor.execute(query("CREATE SCHEMA {s}"))
            cursor.execute(
                query("""CREATE TABLE {s}.entity_evidence_resolution (
                source_id bigint NOT NULL, entity_evidence_id uuid NOT NULL,
                molecular_type_id smallint) PARTITION BY LIST(source_id)""")
            )
            cursor.execute(
                query("""CREATE TABLE {s}.entity_evidence_resolution_default
                PARTITION OF {s}.entity_evidence_resolution DEFAULT""")
            )
            for source in (1, 2):
                cursor.execute(
                    sql.SQL("CREATE TABLE {}.{} PARTITION OF {}.{} FOR VALUES IN ({})").format(
                        sql.Identifier(schema),
                        sql.Identifier(f"{TABLE}_source_{source}"),
                        sql.Identifier(schema),
                        sql.Identifier(TABLE),
                        sql.Literal(source),
                    )
                )
            cursor.executemany(
                query("INSERT INTO {s}.entity_evidence_resolution VALUES(%s,%s,%s)").as_string(
                    connection
                ),
                [
                    (source, str(uuid.UUID(int=ordinal)), molecular)
                    for ordinal, (source, molecular) in enumerate(
                        ((1, 2), (1, None), (2, 2), (99, 7)), 1
                    )
                ],
            )
            cursor.execute(
                query(
                    "CREATE TABLE {s}.parquet_release (singleton boolean PRIMARY KEY, constraint_plan jsonb)"
                )
            )
            cursor.execute(query("INSERT INTO {s}.parquet_release VALUES (true,'{{}}')"))
            cursor.execute(
                query(
                    "CREATE INDEX entity_evidence_resolution_molecular_type_idx "
                    "ON {s}.entity_evidence_resolution(molecular_type_id)"
                )
            )
        connection.commit()
        try:
            yield connection, schema
        finally:
            connection.rollback()
            with connection.cursor() as cursor:
                cursor.execute(query("DROP SCHEMA {s} CASCADE"))
            connection.commit()


def _state(connection, schema):
    with connection.cursor() as cursor:
        cursor.execute(
            """SELECT table_class.relname,index_class.relname,index_class.relkind,
            pg_get_indexdef(index_info.indexrelid),index_info.indisvalid,index_info.indisready,
            parent_class.relname
            FROM pg_index index_info JOIN pg_class table_class ON table_class.oid=index_info.indrelid
            JOIN pg_namespace namespace ON namespace.oid=table_class.relnamespace
            JOIN pg_class index_class ON index_class.oid=index_info.indexrelid
            LEFT JOIN pg_inherits attachment ON attachment.inhrelid=index_class.oid
            LEFT JOIN pg_class parent_class ON parent_class.oid=attachment.inhparent
            WHERE namespace.nspname=%s AND table_class.relname=ANY(%s)
            ORDER BY table_class.relname,index_class.relname""",
            (schema, [TABLE, TABLE + "_default", TABLE + "_source_1", TABLE + "_source_2"]),
        )
        return Counter(
            tuple(
                value.replace(schema + ".", "<schema>.") if isinstance(value, str) else value
                for value in row
            )
            for row in cursor.fetchall()
        )


def _oids(connection, schema):
    with connection.cursor() as cursor:
        cursor.execute(
            """SELECT name.relname,name.oid FROM pg_class name
            JOIN pg_namespace namespace ON namespace.oid=name.relnamespace
            WHERE namespace.nspname=%s AND name.relkind IN ('i','I') ORDER BY name.relname""",
            (schema,),
        )
        return cursor.fetchall()


def _rows(connection, schema):
    with connection.cursor() as cursor:
        cursor.execute(
            sql.SQL("SELECT * FROM {}.{}").format(sql.Identifier(schema), sql.Identifier(TABLE))
        )
        return Counter(cursor.fetchall())


def _assert_complete_valid(state):
    assert len(state) == 4
    assert {row[0] for row in state} == {
        TABLE,
        TABLE + "_default",
        TABLE + "_source_1",
        TABLE + "_source_2",
    }
    assert all(row[4] is True and row[5] is True for row in state)
    assert {row[6] for row in state if row[0] != TABLE} == {INDEX}
    assert [row for row in state if row[0] == TABLE][0][6] is None


def test_captured_partition_parent_replay_restores_exact_children_and_validity(index_database):
    connection, schema = index_database
    original, original_rows = _state(connection, schema), _rows(connection, schema)
    _assert_complete_valid(original)
    plan = aligned_loader._defer_constraints(connection, schema, [TABLE])
    assert len(plan["indexes"]) == 1
    assert plan["indexes"][0][:2] == (TABLE, INDEX)
    # PostgreSQL's introspection serializes even a VALID partitioned parent
    # using ON ONLY. Replaying that literally creates an invalid empty parent.
    assert " ON ONLY " in plan["indexes"][0][2]
    assert not _state(connection, schema)
    aligned_loader._restore_constraints(connection, schema, plan)
    connection.commit()
    restored = _state(connection, schema)
    _assert_complete_valid(restored)
    assert restored == original
    assert _rows(connection, schema) == original_rows


def test_molecular_index_helper_is_idempotent_for_valid_parent_and_children(index_database):
    connection, schema = index_database
    original = _state(connection, schema)
    original_oids, original_rows = _oids(connection, schema), _rows(connection, schema)
    for _ in range(2):
        aligned_loader._ensure_molecular_type_index(connection, schema)
        current = _state(connection, schema)
        _assert_complete_valid(current)
        assert current == original
        assert _oids(connection, schema) == original_oids
        assert _rows(connection, schema) == original_rows


@pytest.mark.parametrize("partially_attached", [False, True])
def test_molecular_index_helper_repairs_real_invalid_parent_without_row_loss(
    index_database, partially_attached
):
    connection, schema = index_database
    original, original_rows = _state(connection, schema), _rows(connection, schema)
    with connection.cursor() as cursor:
        cursor.execute(
            sql.SQL("DROP INDEX {}.{} CASCADE").format(
                sql.Identifier(schema), sql.Identifier(INDEX)
            )
        )
        cursor.execute(
            sql.SQL("CREATE INDEX {} ON ONLY {}.{} (molecular_type_id)").format(
                sql.Identifier(INDEX), sql.Identifier(schema), sql.Identifier(TABLE)
            )
        )
        if partially_attached:
            # Preserve one already-attached child and reuse another equivalent
            # valid index that exists but is not yet attached to the parent.
            by_table = {row[0]: row[1] for row in original}
            for source in (1, 2):
                table = TABLE + f"_source_{source}"
                cursor.execute(
                    sql.SQL("CREATE INDEX {} ON {}.{} (molecular_type_id)").format(
                        sql.Identifier(by_table[table]),
                        sql.Identifier(schema),
                        sql.Identifier(table),
                    )
                )
            cursor.execute(
                sql.SQL("ALTER INDEX {}.{} ATTACH PARTITION {}.{}").format(
                    sql.Identifier(schema),
                    sql.Identifier(INDEX),
                    sql.Identifier(schema),
                    sql.Identifier(by_table[TABLE + "_source_1"]),
                )
            )
    connection.commit()
    broken = _state(connection, schema)
    assert len(broken) == (3 if partially_attached else 1)
    parent = next(row for row in broken if row[0] == TABLE)
    assert parent[0] == TABLE and parent[1] == INDEX
    assert parent[4] is False and parent[5] is True
    assert parent[6] is None
    retained_oids = dict(_oids(connection, schema))
    aligned_loader._ensure_molecular_type_index(connection, schema)
    repaired = _state(connection, schema)
    _assert_complete_valid(repaired)
    assert repaired == original
    assert _rows(connection, schema) == original_rows
    assert all(
        dict(_oids(connection, schema))[name] == value for name, value in retained_oids.items()
    )
    repaired_oids = _oids(connection, schema)
    aligned_loader._ensure_molecular_type_index(connection, schema)
    assert _oids(connection, schema) == repaired_oids
