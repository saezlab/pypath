"""The schema the loader creates holds no table that nothing fills or reads."""

from contextlib import closing
import uuid

import psycopg2
from psycopg2 import sql
import pytest

from omnipath_postgres.relational.db.schema import ensure_schema

DROPPED = (
    "resolver_chemical_identifier_lookup",
    "resolver_chemical_identifier_lookup_ambiguous",
    "resolver_protein_identifier_lookup",
    "resolver_protein_identifier_lookup_ambiguous",
    "gene_protein_representative",
    "state",
    "state_component",
    "evidence_state",
)


@pytest.fixture(scope="module")
def tables(postgres_dsn):
    schema = "schema_" + uuid.uuid4().hex
    with closing(psycopg2.connect(postgres_dsn)) as connection:
        try:
            ensure_schema(connection, schema=schema, indexes=False, progress=False)
            with connection.cursor() as cur:
                cur.execute(
                    "SELECT c.relname FROM pg_class c JOIN pg_namespace n ON n.oid = c.relnamespace "
                    "WHERE n.nspname = %s AND c.relkind IN ('r', 'p', 'v', 'm')",
                    [schema],
                )
                yield {row[0] for row in cur.fetchall()}
        finally:
            connection.rollback()
            with connection.cursor() as cur:
                cur.execute(sql.SQL("DROP SCHEMA IF EXISTS {} CASCADE").format(sql.Identifier(schema)))
            connection.commit()


@pytest.mark.parametrize("table", DROPPED)
def test_a_dead_table_is_not_created(tables, table):
    assert table not in tables


def test_the_record_of_resolution_conflicts_stays(tables):
    assert "resolution_conflict" in tables


def test_the_tables_the_loader_fills_are_still_there(tables):
    assert {"entity", "relation_evidence", "identifier_evidence", "vocab_relation_predicate"} <= tables
