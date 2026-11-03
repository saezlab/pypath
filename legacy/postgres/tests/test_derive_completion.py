"""Integration tests for derive completion + count/overlap tables.

Asserts the post-`derive` state of a built OmniPath database:
- the derive ran end-to-end (resources, facet bitmaps, and the ontology
  entity-relations tables populated),
- the new derived tables ``entity_source_count`` and ``resource_overlap_summary``
  exist, are populated, and answer their queries quickly,
- no *unattached* ``*_staging`` leftover tables remain (the orphan sweep ran;
  attached ``*_staging`` partitions are live data and are kept).

Run against a built instance, e.g. on beauty::

    DATABASE_URL=postgresql://omnipath:omnipath@localhost:5404/omnipath \
        uv run pytest tests/test_derive_completion.py -v

Skipped when DATABASE_URL is not set (no DB to test against).
"""

from __future__ import annotations

import os
import time

import pytest

DATABASE_URL = os.environ.get('DATABASE_URL')
SCHEMA = os.environ.get('OMNIPATH_PG_SCHEMA', 'public')

pytestmark = pytest.mark.skipif(
    not DATABASE_URL,
    reason='DATABASE_URL not set; derive-completion test needs a built database',
)


@pytest.fixture(scope='module')
def conn():
    import psycopg2

    connection = psycopg2.connect(DATABASE_URL)
    try:
        yield connection
    finally:
        connection.close()


def _scalar(conn, query: str):
    with conn.cursor() as cur:
        cur.execute(query)
        return cur.fetchone()[0]


def _table_exists(conn, name: str) -> bool:
    with conn.cursor() as cur:
        cur.execute('SELECT to_regclass(%s)', [f'{SCHEMA}.{name}'])
        return cur.fetchone()[0] is not None


def test_derive_core_tables_populated(conn):
    """resources, facet bitmaps and the ontology tables are populated."""
    assert _scalar(conn, f'SELECT count(*) FROM {SCHEMA}.resources') > 0
    assert (
        _scalar(
            conn,
            f"SELECT count(*) FROM {SCHEMA}.facet_entity_bitmap "
            f"WHERE facet_name = 'source'",
        )
        > 0
    )
    assert (
        _scalar(
            conn,
            f"SELECT count(*) FROM {SCHEMA}.facet_relation_bitmap "
            f"WHERE facet_name = 'source'",
        )
        > 0
    )
    # Ontology now lives in the entity-relations schema (the legacy
    # ``ontology_terms`` table + OBO artifacts + ontograph are deprecated):
    # terms are entities (``Cv Term:OM:0012``), the hierarchy is
    # ``entity_ontology_relation`` and ``entity_ontology_term`` is the
    # denormalised derived table the API reads.
    assert _scalar(conn, f'SELECT count(*) FROM {SCHEMA}.entity_ontology_term') > 0
    assert (
        _scalar(conn, f'SELECT count(*) FROM {SCHEMA}.entity_ontology_relation')
        > 0
    )


def test_entity_source_count(conn):
    """entity_source_count exists, is populated, and the coverage profile is fast."""
    assert _table_exists(conn, 'entity_source_count')
    assert _scalar(conn, f'SELECT count(*) FROM {SCHEMA}.entity_source_count') > 0
    # source_count is positive and source_list length matches source_count.
    bad = _scalar(
        conn,
        f"SELECT count(*) FROM {SCHEMA}.entity_source_count "
        f"WHERE source_count < 1 "
        f"OR source_count <> cardinality(source_list)",
    )
    assert bad == 0

    # Coverage-profile query (Panel B) returns quickly from the derived table.
    started = time.perf_counter()
    with conn.cursor() as cur:
        cur.execute(
            f'SELECT source_count AS n_resources, count(*) AS n_entities '
            f'FROM {SCHEMA}.entity_source_count GROUP BY source_count '
            f'ORDER BY source_count'
        )
        rows = cur.fetchall()
    elapsed = time.perf_counter() - started
    assert rows  # at least one (n_resources, n_entities) bucket
    assert elapsed < 1.0, f'coverage profile took {elapsed:.3f}s (>1s)'


def test_resource_overlap_summary(conn):
    """resource_overlap_summary exists and holds entity + relation overlaps."""
    assert _table_exists(conn, 'resource_overlap_summary')
    kinds = {
        row
        for (row,) in _rows(
            conn,
            f'SELECT DISTINCT content_kind FROM {SCHEMA}.resource_overlap_summary',
        )
    }
    assert 'entity' in kinds
    # bounded: at most N*N per content kind
    n_sources = _scalar(conn, f'SELECT count(*) FROM {SCHEMA}.data_source')
    n_overlap = _scalar(
        conn, f'SELECT count(*) FROM {SCHEMA}.resource_overlap_summary'
    )
    assert n_overlap <= n_sources * n_sources * 3
    # every overlap is positive and each unordered pair stored once
    assert (
        _scalar(
            conn,
            f'SELECT count(*) FROM {SCHEMA}.resource_overlap_summary '
            f'WHERE overlap < 1 OR source_a_id >= source_b_id',
        )
        == 0
    )

    started = time.perf_counter()
    list(
        _rows(
            conn,
            f"SELECT source_a_id, source_b_id, overlap "
            f"FROM {SCHEMA}.resource_overlap_summary "
            f"WHERE content_kind = 'entity' ORDER BY overlap DESC",
        )
    )
    elapsed = time.perf_counter() - started
    assert elapsed < 1.0, f'resource overlap query took {elapsed:.3f}s (>1s)'


def test_no_unattached_staging_tables(conn):
    """No unattached *_staging orphans remain after the sweep.

    The ``*_source_<N>_staging`` tables are the per-source LIST partitions of the
    evidence tables (the loader CREATEs, COPYs, then ``ATTACH PARTITION`` keeping
    the name) and are live data that must be kept. Only a staging table created
    but never attached (an aborted load) is an orphan; the sweep drops those.
    """
    leftover = _scalar(
        conn,
        """
        SELECT count(*)
        FROM pg_class c
        JOIN pg_namespace n ON n.oid = c.relnamespace
        WHERE n.nspname = '%s'
          AND c.relkind = 'r'
          AND c.relname ~ '_source_[0-9]+_staging$'
          AND NOT EXISTS (
            SELECT 1 FROM pg_inherits i WHERE i.inhrelid = c.oid
          )
        """
        % SCHEMA,
    )
    assert leftover == 0, f'{leftover} unattached *_staging tables remain'


def _rows(conn, query: str):
    with conn.cursor() as cur:
        cur.execute(query)
        return cur.fetchall()
