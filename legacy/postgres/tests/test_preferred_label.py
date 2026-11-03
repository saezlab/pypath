"""spec 011 T065 -- every chemical entity carries exactly one preferred label.

The label cascade (``chemical_labels.py``) gained a new pass (T074): the
shortest name ``entity_name`` (T072/T073) marks ``is_preferred`` -- the
minting resource's own recommended name, additive on top of the pre-existing
name cascade. Asserts the post-``derive`` state: ``entity_name`` exists and
is populated, every chemical entity still carries exactly one non-empty
label (the new pass must never leave one unlabeled or double-labeled), and
an entity labeled through the new pass actually matches a real
``is_preferred`` row.

Run against a built instance, e.g. on beauty::

    DATABASE_URL=postgresql://omnipath:omnipath@localhost:55432/omnipath \
        uv run pytest tests/test_preferred_label.py -v

Skipped when DATABASE_URL is not set (no DB to test against).
"""

from __future__ import annotations

import os

import pytest

DATABASE_URL = os.environ.get('DATABASE_URL')
SCHEMA = os.environ.get('OMNIPATH_PG_SCHEMA', 'public')

CHEMICAL_ENTITY_TYPE = 'Chemical:OM:0037'

pytestmark = pytest.mark.skipif(
    not DATABASE_URL,
    reason='DATABASE_URL not set; preferred-label test needs a built database',
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


def _chem_type_id(conn):
    return _scalar(
        conn,
        f'SELECT entity_type_id FROM {SCHEMA}.vocab_entity_type WHERE name = %s',
        [CHEMICAL_ENTITY_TYPE],
    )


def test_entity_name_table_is_populated(conn):
    count = _scalar(conn, f'SELECT count(*) FROM {SCHEMA}.entity_name')
    assert count and count > 0, 'entity_name is empty -- did derive run?'


def test_entity_name_is_preferred_is_never_null(conn):
    # The column is NOT NULL in the schema; a stale pre-migration table
    # would still surface here as a real gap.
    bad = _scalar(
        conn,
        f'SELECT count(*) FROM {SCHEMA}.entity_name WHERE is_preferred IS NULL',
    )
    assert bad == 0


def test_every_chemical_entity_has_exactly_one_label(conn):
    chem = _chem_type_id(conn)
    without_label = _scalar(
        conn,
        f"""
        SELECT count(*) FROM {SCHEMA}.entity
        WHERE entity_type_id = %s
          AND (label IS NULL OR label = '' OR label_rule IS NULL)
        """,
        [chem],
    )
    assert without_label == 0


def test_preferred_name_pass_matches_a_real_preferred_row(conn):
    chem = _chem_type_id(conn)
    mismatched = _scalar(
        conn,
        f"""
        SELECT count(*) FROM {SCHEMA}.entity e
        WHERE e.entity_type_id = %s
          AND e.label_rule = 'chemical_preferred_name'
          AND NOT EXISTS (
            SELECT 1 FROM {SCHEMA}.entity_name en
            WHERE en.entity_id = e.entity_id
              AND en.is_preferred
              AND en.name = e.label
          )
        """,
        [chem],
    )
    assert mismatched == 0
