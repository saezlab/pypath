"""The step accepts a session the caller has already used.

The derive hands MetSigDB a long-lived connection that has already created
temp tables, and Postgres refuses `SET temp_buffers` once a session has touched
one. The refusal arrives as `invalid_parameter_value` and aborts the
transaction, so the step died on its first statement and the derive logged the
failure and carried on. Under `make`, which sets a quiet logging config, that
error was never printed: the build reported success and published nothing.

A fresh connection never reproduces it, which is why the standalone build, the
quickstart and 32 build-side tests all passed while every derive silently
skipped the step.

    DATABASE_URL=postgresql://omnipath:omnipath@localhost:55435/omnipath \
        uv run --with pytest --with psycopg2-binary pytest \
        tests/test_metsigdb_reused_session.py -v
"""

from __future__ import annotations

import os

import pytest

DATABASE_URL = os.environ.get('DATABASE_URL')

pytestmark = pytest.mark.skipif(
    not DATABASE_URL,
    reason='DATABASE_URL not set; this exercises a real session',
)


@pytest.fixture
def used_conn():
    """A connection that has touched a temp table, as the derive's has."""
    import psycopg2

    connection = psycopg2.connect(DATABASE_URL)
    with connection.cursor() as cur:
        cur.execute('CREATE TEMP TABLE metsigdb_probe (x int) ON COMMIT DROP')
        cur.execute('INSERT INTO metsigdb_probe VALUES (1)')
    try:
        yield connection
    finally:
        connection.rollback()
        connection.close()


def test_widening_temp_buffers_survives_a_used_session(used_conn):
    """The refusal is tolerated, and it does not poison the transaction."""
    from omnipath_build.metsigdb.build import _widen_temp_buffers

    with used_conn.cursor() as cur:
        _widen_temp_buffers(cur)
        # The real regression: before the fix the aborted transaction made
        # every following statement raise InFailedSqlTransaction, so the step
        # failed on whatever it tried next rather than on the SET.
        cur.execute('SELECT 1')
        assert cur.fetchone()[0] == 1


def test_widening_temp_buffers_still_applies_on_a_fresh_session():
    """On a session that has touched no temp table, the SET must still take."""
    import psycopg2

    from omnipath_build.metsigdb.build import TEMP_BUFFERS, _widen_temp_buffers

    connection = psycopg2.connect(DATABASE_URL)
    try:
        with connection.cursor() as cur:
            _widen_temp_buffers(cur)
            cur.execute('SHOW temp_buffers')
            applied = cur.fetchone()[0]
        assert applied.upper().replace(' ', '') == TEMP_BUFFERS.upper()
    finally:
        connection.rollback()
        connection.close()
