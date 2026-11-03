"""Integration tests for the build manifest.

Run against a built instance after `derive`, e.g. on beauty::

    MAX_RECORDS=1000 \
    DATABASE_URL=postgresql://omnipath:omnipath@localhost:5404/omnipath \
        uv run --with pytest pytest tests/test_build_manifest.py -v

The instance is expected to be a capped build (MAX_RECORDS set), so
`partial_build` is asserted true. Skipped when DATABASE_URL is not set.
"""

from __future__ import annotations

import os
import re

import pytest

DATABASE_URL = os.environ.get('DATABASE_URL')
SCHEMA = os.environ.get('OMNIPATH_PG_SCHEMA', 'public')

pytestmark = pytest.mark.skipif(
    not DATABASE_URL,
    reason='DATABASE_URL not set; build-manifest test needs a built database',
)


@pytest.fixture(scope='module')
def conn():
    import psycopg2

    connection = psycopg2.connect(DATABASE_URL)
    try:
        yield connection
    finally:
        connection.close()


def _row(conn, query, params=None):
    with conn.cursor() as cur:
        cur.execute(query, params)
        return cur.fetchone()


def test_exactly_one_manifest_row(conn):
    """A fresh database carries exactly one current build_manifest row."""
    count = _row(conn, f'SELECT count(*) FROM {SCHEMA}.build_manifest')[0]
    assert count == 1


def test_manifest_fields_well_formed(conn):
    """All required fields are present and well-formed."""
    build_id, built_at, package_commits, resources, partial_build = _row(
        conn,
        f"""
        SELECT build_id, built_at, package_commits, resources, partial_build
        FROM {SCHEMA}.build_manifest
        """,
    )
    assert re.fullmatch(r'[0-9a-f]{12}', build_id)
    assert built_at is not None
    assert isinstance(package_commits, dict)
    assert {'omnipath_build', 'omnipath_resources'} <= set(package_commits)
    for pkg in package_commits.values():
        assert 'commit' in pkg and 'dirty' in pkg
    assert isinstance(resources, list) and resources
    for item in resources:
        assert 'name' in item and 'record_count' in item
        assert isinstance(item['record_count'], int)
    assert isinstance(partial_build, bool)


def test_partial_build_flagged_on_capped_build(conn):
    """A MAX_RECORDS-capped build sets partial_build true."""
    partial_build = _row(
        conn, f'SELECT partial_build FROM {SCHEMA}.build_manifest'
    )[0]
    assert partial_build is True


def test_build_id_reproducible(conn):
    """build_id is the 12-hex content hash of {package_commits, resources}."""
    import hashlib
    import json

    package_commits, resources, build_id = _row(
        conn,
        f"""
        SELECT package_commits, resources, build_id
        FROM {SCHEMA}.build_manifest
        """,
    )
    payload = {'package_commits': package_commits, 'resources': resources}
    expected = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(',', ':')).encode()
    ).hexdigest()[:12]
    assert build_id == expected
