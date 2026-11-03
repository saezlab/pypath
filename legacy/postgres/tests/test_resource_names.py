"""Build-side tests for the resource 3-name model.

Run against a built instance after `derive`, e.g. on beauty::

    DATABASE_URL=postgresql://omnipath:omnipath@localhost:5404/omnipath \
        uv run --with pytest pytest tests/test_resource_names.py -v

Skipped when DATABASE_URL is not set.
"""

from __future__ import annotations

import os
import re

import pytest

DATABASE_URL = os.environ.get('DATABASE_URL')
SCHEMA = os.environ.get('OMNIPATH_PG_SCHEMA', 'public')

pytestmark = pytest.mark.skipif(
    not DATABASE_URL, reason='DATABASE_URL not set; resource-name test needs a build'
)

_SLUG_RE = re.compile(r'^[a-z0-9]+$')


@pytest.fixture(scope='module')
def conn():
    import psycopg2

    connection = psycopg2.connect(DATABASE_URL)
    try:
        yield connection
    finally:
        connection.close()


def _rows(conn, query):
    with conn.cursor() as cur:
        cur.execute(query)
        return cur.fetchall()


def test_resources_have_short_and_full(conn):
    """Every resource carries resource_short + resource_full (3-name model)."""
    missing = _rows(
        conn,
        f'SELECT resource_id FROM {SCHEMA}.resources '
        f'WHERE resource_short IS NULL OR resource_full IS NULL',
    )
    assert missing == []


def test_json_backed_resources_use_authoritative_names(conn):
    """Resources present in resources.json resolve to their authoritative names."""
    names = dict(
        (rid, (short, full))
        for rid, short, full in _rows(
            conn,
            f'SELECT resource_id, resource_short, resource_full '
            f'FROM {SCHEMA}.resources',
        )
    )
    if 'signor' in names:
        assert names['signor'] == ('SIGNOR', 'Signaling Network Open Resource')
    if 'chembl' in names:
        assert names['chembl'][0] == 'ChEMBL'
    if 'cellphonedb' in names:
        assert names['cellphonedb'][0] == 'CellPhoneDB'


def test_no_underscore_in_names(conn):
    """short / full contain no underscore (the reserved primary_secondary char).

    Spaces in inputs_v2-only resources' short names (not in the legacy
    resources.json) are flagged by the build validator and curated per-module in
    the resource-metadata audit; here we only assert the underscore invariant
    holds.
    """
    for rid, short, full in _rows(
        conn,
        f'SELECT resource_id, resource_short, resource_full FROM {SCHEMA}.resources '
        f'WHERE resource_short IS NOT NULL',
    ):
        assert '_' not in (short or ''), f'{rid} short has _'
        assert '_' not in (full or ''), f'{rid} full has _'
