"""spec 011 T125 -- the manifest records chemical resolution coverage, and
two builds compare from the manifests alone.

Run against a built instance, e.g. on beauty::

    DATABASE_URL=postgresql://omnipath:omnipath@localhost:55432/omnipath \
        uv run pytest tests/test_coverage_manifest.py -v

Skipped when DATABASE_URL is not set (no DB to test against).
"""

from __future__ import annotations

import os

import pytest

DATABASE_URL = os.environ.get('DATABASE_URL')
SCHEMA = os.environ.get('OMNIPATH_PG_SCHEMA', 'public')

pytestmark = pytest.mark.skipif(
    not DATABASE_URL,
    reason='DATABASE_URL not set; the coverage manifest needs a built database',
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


def test_chemical_resolution_coverage_table_exists(conn):
    present = _scalar(
        conn,
        f"SELECT to_regclass('{SCHEMA}.chemical_resolution_coverage')",
    )
    assert present is not None


def test_coverage_is_populated(conn):
    from omnipath_build.db.derived_tables import (
        _populate_chemical_resolution_coverage,
    )

    with conn.cursor() as cur:
        rows = _populate_chemical_resolution_coverage(cur, SCHEMA)
    conn.commit()
    assert rows > 0

    count = _scalar(
        conn, f'SELECT count(*) FROM {SCHEMA}.chemical_resolution_coverage',
    )
    assert count == rows


def test_coverage_rows_are_internally_consistent(conn):
    """entities <= mentions, and each outcome bucket <= mentions."""
    rows = None
    with conn.cursor() as cur:
        cur.execute(
            f"""
            SELECT mentions, entities, reached_structure, reached_name,
                   unresolved, conflicted
            FROM {SCHEMA}.chemical_resolution_coverage
            """
        )
        rows = cur.fetchall()
    assert rows, 'no coverage rows -- did the populate step run?'
    for mentions, entities, structure, name, unresolved, conflicted in rows:
        assert entities <= mentions
        assert structure <= mentions
        assert name <= mentions
        assert unresolved <= mentions
        assert conflicted <= mentions


def test_manifest_records_chemical_resolution_coverage(conn):
    coverage = _scalar(
        conn,
        f'SELECT chemical_resolution_coverage FROM {SCHEMA}.build_manifest '
        'LIMIT 1',
    )
    assert coverage is not None, (
        'build_manifest.chemical_resolution_coverage is NULL -- '
        'did resources.py record it?'
    )


def test_two_builds_compare_from_the_manifest_alone(conn):
    """The manifest's own JSON is enough to diff coverage across builds --
    no re-running analysis, no re-querying content tables.
    """
    manifest = _scalar(
        conn,
        f'SELECT chemical_resolution_coverage FROM {SCHEMA}.build_manifest '
        'LIMIT 1',
    )
    assert isinstance(manifest, list)
    assert manifest, 'manifest coverage list is empty'
    row = manifest[0]
    for key in (
        'source_id', 'identifier_type_id', 'role', 'mentions', 'entities',
        'reached_structure', 'reached_name', 'unresolved', 'conflicted',
    ):
        assert key in row
