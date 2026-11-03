"""Integration test for the MAX_RECORDS capped full cycle.

Two layers:

1. A light invariant (runs whenever DATABASE_URL points at a MAX_RECORDS build):
   the build advertises itself as partial — capped builds are never silently
   treated as authoritative.

2. An opt-in heavy benchmark (only when OMNIPATH_RUN_FULL_BUILD=1): run a full
   capped `make all` cycle and assert it finishes under the 10-minute budget and
   sets partial_build. Heavy + environment-specific, so it is skipped by default.

Run the heavy layer on beauty, e.g.::

    OMNIPATH_RUN_FULL_BUILD=1 OMNIPATH_BUILD_MAX_RECORDS=1000 \
    OMNIPATH_BUILD_DIR=~/instances/dev4/omnipath-build \
    DATABASE_URL=postgresql://omnipath:omnipath@localhost:5404/omnipath \
        uv run --with pytest pytest tests/test_max_records.py -v
"""

from __future__ import annotations

import os
import subprocess
import time

import pytest

DATABASE_URL = os.environ.get('DATABASE_URL')
SCHEMA = os.environ.get('OMNIPATH_PG_SCHEMA', 'public')
MAX_RECORDS = os.environ.get('OMNIPATH_BUILD_MAX_RECORDS', '1000')
BUILD_DIR = os.environ.get('OMNIPATH_BUILD_DIR', '.')
CYCLE_BUDGET_SECONDS = 10 * 60


@pytest.fixture(scope='module')
def conn():
    if not DATABASE_URL:
        pytest.skip('DATABASE_URL not set; MAX_RECORDS test needs a built database')
    import psycopg2

    connection = psycopg2.connect(DATABASE_URL)
    try:
        yield connection
    finally:
        connection.close()


def _scalar(conn, query):
    with conn.cursor() as cur:
        cur.execute(query)
        return cur.fetchone()[0]


def test_capped_build_is_flagged_partial(conn):
    """A MAX_RECORDS build advertises partial_build — never silently authoritative."""
    partial = _scalar(
        conn, f'SELECT partial_build FROM {SCHEMA}.build_manifest'
    )
    assert partial is True


@pytest.mark.skipif(
    os.environ.get('OMNIPATH_RUN_FULL_BUILD') != '1',
    reason='heavy: set OMNIPATH_RUN_FULL_BUILD=1 to run the full capped cycle',
)
def test_capped_full_cycle_under_budget(conn):
    """A full capped `make all` cycle finishes under budget and flags partial."""
    env = {**os.environ, 'MAX_RECORDS': MAX_RECORDS, 'DROP_EXISTING': '1'}
    started = time.perf_counter()
    result = subprocess.run(
        ['make', 'all'],
        cwd=BUILD_DIR,
        env=env,
        timeout=CYCLE_BUDGET_SECONDS,
        capture_output=True,
        text=True,
    )
    elapsed = time.perf_counter() - started
    assert result.returncode == 0, result.stderr[-2000:]
    assert elapsed < CYCLE_BUDGET_SECONDS, f'cycle took {elapsed:.0f}s'
    assert _scalar(conn, f'SELECT partial_build FROM {SCHEMA}.build_manifest') is True
