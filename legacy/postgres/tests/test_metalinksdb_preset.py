"""MetaLinksDB as a preset: what must be gone, and what must not be.

Two halves, and the second is the one that is easy to get wrong.

**Gone**: the fifteen materialized views. A dataset expressed as a preset that
still carries a view of its own has not been converted, it has been duplicated,
and the duplicate is the copy that goes stale.

**Not gone**: the rows the curation excludes. Curation is a *query scope*, not
a deletion. ChEMBL's non-mechanism assertions and every BindingDB row stay in
the interaction record, reachable by any other query, because a resource this
dataset does not want is still a resource the database holds. A build that
enforced the curation by dropping rows would make the dataset cheap and the
database wrong.

The second half has a failure mode the first has not: it passes for a build in
which the resource was never loaded. So each assertion first establishes that
the resource contributes at all, and says so when it does not, rather than
reporting an absent resource as a curation success.

Run against a build database, e.g. on dev4::

    DATABASE_URL=postgresql://omnipath:omnipath@localhost:5404/omnipath \
        uv run --with pytest pytest tests/test_metalinksdb_preset.py -v
"""

from __future__ import annotations

import os

import pytest

DATABASE_URL = os.environ.get('DATABASE_URL')
SCHEMA = os.environ.get('OMNIPATH_PG_SCHEMA', 'public')

pytestmark = pytest.mark.skipif(
    not DATABASE_URL,
    reason='DATABASE_URL not set; the preset conversion needs a built database',
)

PRESET = 'metalinksdb'

# The prefix every bespoke view of this dataset carried. Fifteen of them, and
# the assertion is that none is left rather than that a particular list is.
VIEW_PREFIX = 'metalinksdb_'

# Restricted to its mechanism-of-action assertions inside the dataset, and
# whole outside it.
CURATED = 'chembl'

# Excluded from the dataset, retained in the record.
EXCLUDED = 'bindingdb'

# The annotation term the mechanism restriction selects on.
MECHANISM_TERM = 'Chembl Mechanism:OM:0227'


@pytest.fixture(scope='module')
def conn():
    """An open connection to the built database."""
    psycopg2 = pytest.importorskip('psycopg2')
    connection = psycopg2.connect(DATABASE_URL)
    try:
        yield connection
    finally:
        connection.close()


def _count(conn, statement, params=None):
    """One scalar count.

    Args:
        conn: An open connection.
        statement: A statement returning one number.
        params: Its parameters.

    Returns:
        The number.
    """
    with conn.cursor() as cur:
        cur.execute(statement, params or [])
        return cur.fetchone()[0]


def _source_rows(conn, name):
    """Record rows one resource contributes, whatever the dataset wants.

    Args:
        conn: An open connection.
        name: The resource name.

    Returns:
        The row count.
    """
    return _count(
        conn,
        f"""
        SELECT count(*) FROM {SCHEMA}.interaction_fact_resource f
        JOIN {SCHEMA}.data_source d USING (source_id)
        WHERE d.name = %s
        """,
        [name],
    )


def test_no_bespoke_view_of_this_dataset_survives(conn):
    """The preset replaces the views; it does not stand beside them."""
    with conn.cursor() as cur:
        cur.execute(
            'SELECT schemaname, matviewname FROM pg_matviews '
            'WHERE matviewname LIKE %s ORDER BY 2',
            [f'{VIEW_PREFIX}%'],
        )
        remaining = [f'{schema}.{name}' for schema, name in cur.fetchall()]
    assert not remaining, (
        f'{len(remaining)} bespoke views survive the conversion: {remaining}. '
        f'A dataset served twice is served from whichever copy the caller '
        f'happens to reach, and only one of them is refreshed'
    )


def test_the_dataset_declares_no_relation_of_its_own(conn):
    """Its registry row is a preset row: no schema, no combined relation."""
    with conn.cursor() as cur:
        cur.execute(
            f'SELECT schema_name, combined_relation FROM {SCHEMA}.network_registry '
            f'WHERE name = %s',
            [PRESET],
        )
        row = cur.fetchone()
    assert row is not None, f'{PRESET} is not registered'
    schema_name, combined_relation = row
    assert schema_name is None and combined_relation is None, (
        f'{PRESET} still names {schema_name}.{combined_relation}; a preset '
        f'that names a relation is a matview network wearing preset columns'
    )


def test_the_excluded_resource_keeps_its_rows(conn):
    """Excluded from the dataset is not deleted from the database."""
    rows = _source_rows(conn, EXCLUDED)
    assert rows, (
        f'{EXCLUDED} contributes no record row at all. Either it was never '
        f'loaded — in which case this test proves nothing about the curation — '
        f'or the curation was enforced by deleting it, which is the failure '
        f'this test exists to catch'
    )


def test_the_curated_resource_keeps_its_uncurated_rows(conn):
    """ChEMBL's non-mechanism assertions stay in the record."""
    total = _source_rows(conn, CURATED)
    assert total, (
        f'{CURATED} contributes no record row; the curation cannot be shown '
        f'to be a scope rather than a filter on the build'
    )
    mechanism = _count(
        conn,
        f"""
        SELECT count(*)
        FROM {SCHEMA}.relation_evidence_annotation rea
        JOIN {SCHEMA}.annotation a USING (annotation_key)
        JOIN {SCHEMA}.data_source d ON d.source_id = rea.source_id
        WHERE a.term = %s AND d.name = %s
        """,
        [MECHANISM_TERM, CURATED],
    )
    assert mechanism, (
        f'no {CURATED} evidence carries {MECHANISM_TERM}; the restriction has '
        f'nothing to select on'
    )
    assert total > mechanism, (
        f'{CURATED} contributes {total} record rows and {mechanism} mechanism '
        f'annotations. With no uncurated remainder the curation is '
        f'indistinguishable from the resource, and this test cannot tell '
        f'a scope from a deletion'
    )
