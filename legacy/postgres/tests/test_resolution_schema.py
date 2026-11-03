"""Schema-addition assertions for the gene-anchored entity model.

Constitution III: assert each new table/column the gene-anchored model introduces
exists after the schema is created. Self-contained — it builds the schema into a
throwaway scratch schema (no data build needed), so it runs against any Postgres.

Run::

    DATABASE_URL=postgresql://omnipath:omnipath@localhost:5404/omnipath \
        pytest tests/test_resolution_schema.py
"""

from __future__ import annotations

import os

import pytest

DATABASE_URL = os.environ.get('DATABASE_URL')
SCRATCH = os.environ.get('OMNIPATH_TEST_SCRATCH_SCHEMA', 'ga_schema_test')


@pytest.fixture(scope='module')
def scratch_schema():
    if not DATABASE_URL:
        pytest.skip('DATABASE_URL not set; schema test needs a Postgres')
    import psycopg2

    from omnipath_build.db import schema as build_schema

    conn = psycopg2.connect(DATABASE_URL)
    try:
        build_schema.ensure_schema(conn, schema=SCRATCH, drop_existing=True)
        conn.commit()
        yield conn
    finally:
        with conn.cursor() as cur:
            cur.execute(f'DROP SCHEMA IF EXISTS {SCRATCH} CASCADE')
        conn.commit()
        conn.close()


def _has_table(conn, name: str) -> bool:
    with conn.cursor() as cur:
        cur.execute(f"SELECT to_regclass('{SCRATCH}.{name}')")
        return cur.fetchone()[0] is not None


def _has_column(conn, table: str, column: str) -> bool:
    with conn.cursor() as cur:
        cur.execute(
            'SELECT 1 FROM information_schema.columns '
            'WHERE table_schema=%s AND table_name=%s AND column_name=%s',
            (SCRATCH, table, column),
        )
        return cur.fetchone() is not None


@pytest.mark.parametrize(
    'table',
    [
        'vocab_molecular_type',
        'gene_protein_representative',
        'state',
        'state_component',
        'evidence_state',
        'evidence_state_default',
    ],
)
def test_gene_anchored_tables_exist(scratch_schema, table):
    assert _has_table(scratch_schema, table), f'missing table {table}'


@pytest.mark.parametrize(
    'column',
    ['label', 'label_rule', 'resolution_mechanism', 'resolution_detail'],
)
def test_entity_label_and_mechanism_columns(scratch_schema, column):
    assert _has_column(scratch_schema, 'entity', column), (
        f'entity missing column {column}'
    )


def test_resolution_has_molecular_type(scratch_schema):
    assert _has_column(
        scratch_schema, 'entity_evidence_resolution', 'molecular_type_id'
    )


def test_molecular_type_vocab_seeded(scratch_schema):
    with scratch_schema.cursor() as cur:
        cur.execute(f'SELECT name FROM {SCRATCH}.vocab_molecular_type')
        names = {r[0] for r in cur.fetchall()}
    assert {'gene', 'protein', 'mrna', 'mirna', 'lncrna'} <= names
