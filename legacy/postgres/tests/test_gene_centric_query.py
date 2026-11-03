"""The gene-centric output reads no state tables.

The common gene-level query path MUST NOT touch ``state`` / ``state_component`` /
``evidence_state``. These tests assert that from the *query plan*:
- ``gene_output`` (gene + representative UniProt) plans without any state table;
- filtering "records that are proteins" via the indexed ``molecular_type_id`` on
  ``entity_evidence_resolution`` plans without any state table;
- and the gene-centric row shape is present (gene / NCBI Gene id / UniProt), with
  the EGFR benchmark reproducing the classic gene+UniProt layout.

Run against a built instance, e.g. on beauty::

    DATABASE_URL=postgresql://omnipath:omnipath@localhost:55432/omnipath \
        uv run pytest tests/test_gene_centric_query.py -v

Skipped when DATABASE_URL is not set (no DB to test against).
"""

from __future__ import annotations

import json
import os

import pytest

DATABASE_URL = os.environ.get('DATABASE_URL')
SCHEMA = os.environ.get('OMNIPATH_PG_SCHEMA', 'public')

STATE_TABLES = {'state', 'state_component', 'evidence_state'}

pytestmark = pytest.mark.skipif(
    not DATABASE_URL,
    reason='DATABASE_URL not set; gene-centric query test needs a built database',
)


@pytest.fixture(scope='module')
def conn():
    import psycopg2

    connection = psycopg2.connect(DATABASE_URL)
    try:
        yield connection
    finally:
        connection.close()


def _scalar(conn, query: str, params=None):
    with conn.cursor() as cur:
        cur.execute(query, params or [])
        row = cur.fetchone()
        return row[0] if row else None


def _relation_exists(conn, name: str) -> bool:
    return _scalar(conn, 'SELECT to_regclass(%s)', [f'{SCHEMA}.{name}']) is not None


def _plan_relations(conn, query: str) -> set[str]:
    """Return the set of base relation names referenced by a query's plan."""
    with conn.cursor() as cur:
        cur.execute(f'EXPLAIN (FORMAT JSON) {query}')
        plan = cur.fetchone()[0]
    if isinstance(plan, str):
        plan = json.loads(plan)
    names: set[str] = set()

    def walk(node):
        if isinstance(node, dict):
            rel = node.get('Relation Name')
            if rel:
                names.add(rel)
            for value in node.values():
                walk(value)
        elif isinstance(node, list):
            for item in node:
                walk(item)

    walk(plan)
    return names


@pytest.fixture(scope='module')
def require_gene_output(conn):
    if not _relation_exists(conn, 'gene_output'):
        pytest.skip('gene_output view absent — build predates the view')


def test_gene_output_has_gene_centric_columns(conn, require_gene_output):
    with conn.cursor() as cur:
        cur.execute(f'SELECT * FROM {SCHEMA}.gene_output LIMIT 0')
        columns = {desc[0] for desc in cur.description}
    for required in ('gene', 'ncbi_gene_id', 'uniprot'):
        assert required in columns, f'gene_output missing {required!r} column'


def test_gene_output_plan_touches_no_state(conn, require_gene_output):
    """A gene_output read does not reference any state table."""
    relations = _plan_relations(conn, f'SELECT * FROM {SCHEMA}.gene_output')
    leaked = relations & STATE_TABLES
    assert not leaked, f'gene_output plan touches state tables: {leaked}'


def test_protein_filter_plan_touches_no_state(conn):
    """'Records that are proteins' uses molecular_type, not state."""
    if not _relation_exists(conn, 'entity_evidence_resolution'):
        pytest.skip('entity_evidence_resolution absent')
    relations = _plan_relations(
        conn,
        f"""
        SELECT eer.entity_id, count(*)
        FROM {SCHEMA}.entity_evidence_resolution eer
        WHERE eer.molecular_type_id = 2
        GROUP BY eer.entity_id
        """,
    )
    leaked = relations & STATE_TABLES
    assert not leaked, f'protein-filter plan touches state tables: {leaked}'


def test_egfr_gene_centric_row(conn, require_gene_output):
    """EGFR reproduces the classic gene + UniProt layout, with no state join."""
    with conn.cursor() as cur:
        cur.execute(
            f"""
            SELECT gene, ncbi_gene_id, uniprot
            FROM {SCHEMA}.gene_output
            WHERE ncbi_gene_id = '1956'
            """
        )
        row = cur.fetchone()
    if row is None:
        pytest.skip('EGFR gene entity absent in this (capped) build')
    gene, ncbi_gene_id, uniprot = row
    assert ncbi_gene_id == '1956'
    assert uniprot == 'P00533', f'EGFR representative UniProt = {uniprot!r}'
    # gene is the symbol once labels are populated; NCBI id fallback else.
    assert gene in {'EGFR', '1956'}, f'unexpected EGFR gene label {gene!r}'
