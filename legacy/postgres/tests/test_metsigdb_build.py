"""Build-side tests for the MetSigDB membership substrate (cycle 010).

Run against a built instance::

    DATABASE_URL=postgresql://omnipath:omnipath@localhost:55435/omnipath \
        uv run --with pytest pytest tests/test_metsigdb_build.py -v

Skipped when DATABASE_URL is not set.
"""

from __future__ import annotations

import os

import pytest

DATABASE_URL = os.environ.get('DATABASE_URL')

pytestmark = pytest.mark.skipif(
    not DATABASE_URL,
    reason='DATABASE_URL not set; the MetSigDB substrate needs a built database',
)

TABLE = 'metsigdb_membership'

# The shared filter columns the serving layer is allowed to filter on. The
# contract fixes them, and each one carries an index.
FILTER_COLUMNS = (
    'resource',
    'set_type',
    'organism',
    'set_source_id',
    'metabolite_entity_id',
)


@pytest.fixture(scope='module')
def conn():
    import psycopg2

    connection = psycopg2.connect(DATABASE_URL)
    # Autocommit, so one failing assertion query does not poison the
    # connection for every test after it.
    connection.autocommit = True
    try:
        yield connection
    finally:
        connection.close()


def _rows(conn, query, params=None):
    with conn.cursor() as cur:
        cur.execute(query, params)
        return cur.fetchall()


def test_membership_table_exists(conn):
    """The build publishes the substrate under its contracted name."""
    rows = _rows(
        conn,
        'SELECT to_regclass(%s) IS NOT NULL',
        [f'public.{TABLE}'],
    )
    assert rows[0][0], f'{TABLE} is absent; run the MetSigDB build step'


def test_membership_table_columns(conn):
    """Every contract field is a column, and the row carries its build stamp."""
    columns = {
        name
        for (name,) in _rows(
            conn,
            'SELECT column_name FROM information_schema.columns '
            'WHERE table_schema = %s AND table_name = %s',
            ['public', TABLE],
        )
    }
    expected = {
        'metabolite_entity_id',
        'metabolite_label',
        'metabolite_entity_type',
        'inchikey',
        'smiles',
        'hmdb',
        'pubchem',
        'chebi',
        'kegg',
        'resource',
        'set_source_id',
        'set_label',
        'set_type',
        'organism',
        'set_size',
        'set_context',
        'provenance_source',
        'provenance_record',
        'build_id',
    }
    assert expected <= columns, f'missing columns: {sorted(expected - columns)}'
    # `inchi` left the contract: no InChI identifier type exists in the schema.
    assert 'inchi' not in columns


def test_membership_row_identity(conn):
    """Row identity is the primary key, so a rebuild refreshes in place."""
    rows = _rows(
        conn,
        """
        SELECT a.attname
        FROM pg_index i
        JOIN pg_attribute a ON a.attrelid = i.indrelid AND a.attnum = ANY (i.indkey)
        WHERE i.indrelid = %s::regclass AND i.indisprimary
        ORDER BY a.attname
        """,
        [f'public.{TABLE}'],
    )
    assert [name for (name,) in rows] == [
        'metabolite_entity_id',
        'resource',
        'set_source_id',
    ]


def test_filter_columns_are_indexed(conn):
    """The serving layer filters on five columns, and each one has an index."""
    indexed = {
        name
        for (name,) in _rows(
            conn,
            """
            SELECT DISTINCT a.attname
            FROM pg_index i
            JOIN pg_attribute a ON a.attrelid = i.indrelid AND a.attnum = ANY (i.indkey)
            WHERE i.indrelid = %s::regclass
            """,
            [f'public.{TABLE}'],
        )
    }
    missing = set(FILTER_COLUMNS) - indexed
    assert not missing, f'unindexed filter columns: {sorted(missing)}'


def test_set_type_is_constrained(conn):
    """Only the three supported set semantics may be stored."""
    with conn.cursor() as cur:
        cur.execute(
            "SELECT count(*) FROM pg_constraint "
            "WHERE conrelid = %s::regclass AND contype = 'c' "
            "AND pg_get_constraintdef(oid) ILIKE %s",
            [f'public.{TABLE}', '%set_type%'],
        )
        assert cur.fetchone()[0] >= 1


def test_ddl_is_idempotent(conn):
    """Applying the DDL twice changes nothing, so a rebuild is safe."""
    from omnipath_build.metsigdb import ensure_membership_table

    before = _rows(
        conn,
        'SELECT column_name, data_type FROM information_schema.columns '
        'WHERE table_schema = %s AND table_name = %s ORDER BY column_name',
        ['public', TABLE],
    )
    ensure_membership_table(conn)
    ensure_membership_table(conn)
    after = _rows(
        conn,
        'SELECT column_name, data_type FROM information_schema.columns '
        'WHERE table_schema = %s AND table_name = %s ORDER BY column_name',
        ['public', TABLE],
    )
    assert before == after


def _fingerprint(conn):
    """Row count plus a content hash over the row identity and set_size."""
    return _rows(
        conn,
        f"""
        SELECT count(*), md5(string_agg(
          resource || '|' || set_source_id || '|' || metabolite_entity_id::text
          || '|' || set_size::text,
          ',' ORDER BY resource, set_source_id, metabolite_entity_id))
        FROM {TABLE}
        """,
    )[0]


@pytest.fixture(scope='module')
def loaded(conn):
    if not _rows(conn, f'SELECT count(*) FROM {TABLE}')[0][0]:
        pytest.skip('the substrate is empty; run the MetSigDB build first')
    return True


def test_rebuild_is_idempotent(conn, loaded):
    """A second load of one resource changes neither count nor content."""
    from omnipath_build.metsigdb import build_id, load_resource, rule_for

    before = _fingerprint(conn)
    stats = load_resource(conn, rule_for('MACdb'), stamp=build_id(conn))
    after = _fingerprint(conn)
    assert before == after
    assert stats.removed == 0


def test_every_row_carries_the_manifest_build_stamp(conn, loaded):
    stamps = _rows(conn, f'SELECT DISTINCT build_id FROM {TABLE}')
    manifest = _rows(conn, 'SELECT build_id FROM build_manifest')
    assert [s for (s,) in stamps] == [manifest[0][0]]


def test_a_capped_run_holds_only_what_it_loaded(conn, loaded):
    """The cap is for the loop, so the substrate must not keep the full load.

    The resource is restored at the end, because the other tests read a
    complete substrate.
    """
    from omnipath_build.metsigdb import build_id, load_resource, rule_for

    rule = rule_for('MACdb')
    stamp = build_id(conn)
    full = _rows(conn, f"SELECT count(*) FROM {TABLE} WHERE resource = 'MACdb'")[0][0]
    try:
        capped = load_resource(conn, rule, stamp=stamp, max_records=500)
        assert capped.rows == 500
        assert capped.removed == full - 500
        held = _rows(
            conn, f"SELECT count(*) FROM {TABLE} WHERE resource = 'MACdb'"
        )[0][0]
        assert held == 500
    finally:
        restored = load_resource(conn, rule, stamp=stamp)
        assert restored.rows == full


def test_set_size_follows_the_capped_population(conn, loaded):
    """`set_size` counts what the build published, never an upstream total."""
    rows = _rows(
        conn,
        f"""
        SELECT resource, set_source_id, min(set_size), count(*)
        FROM {TABLE} GROUP BY 1, 2
        HAVING min(set_size) <> count(*)
        """,
    )
    assert rows == []


def test_a_build_leaves_the_table_analyzed(conn, loaded):
    """The build replaces every row, so stale statistics misplan every filter."""
    rows = _rows(
        conn,
        "SELECT last_analyze IS NOT NULL OR last_autoanalyze IS NOT NULL "
        "FROM pg_stat_user_tables WHERE relname = %s",
        [TABLE],
    )
    assert rows and rows[0][0]
