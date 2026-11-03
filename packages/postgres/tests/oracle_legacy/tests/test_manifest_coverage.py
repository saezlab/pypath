"""The build manifest records which identifier-resolution tables the build used
and how completely evidence was canonicalised — without letting those volatile
numbers change the build's identity hash.

Run against a built instance after `derive`, e.g. on beauty::

    MAX_RECORDS=1000 \
    DATABASE_URL=postgresql://omnipath:omnipath@localhost:5404/omnipath \
        uv run --with pytest pytest tests/test_manifest_coverage.py -v

Skipped when DATABASE_URL is not set.
"""

from __future__ import annotations

import hashlib
import json
import os

import pytest

DATABASE_URL = os.environ.get('DATABASE_URL')
SCHEMA = os.environ.get('OMNIPATH_PG_SCHEMA', 'public')

pytestmark = pytest.mark.skipif(
    not DATABASE_URL,
    reason='DATABASE_URL not set; manifest-coverage test needs a built database',
)


@pytest.fixture(scope='module')
def conn():
    import psycopg2

    connection = psycopg2.connect(DATABASE_URL)
    try:
        yield connection
    finally:
        connection.close()


def _row(conn, query):
    with conn.cursor() as cur:
        cur.execute(query)
        return cur.fetchone()


def test_manifest_lists_translation_tables(conn):
    translation_tables = _row(
        conn, f'SELECT translation_tables FROM {SCHEMA}.build_manifest'
    )[0]
    assert isinstance(translation_tables, list) and translation_tables
    names = {entry['name'] for entry in translation_tables}
    assert {'resolver_gene', 'resolver_protein', 'resolver_chemical'} <= names


def test_manifest_reports_canonicalization_coverage(conn):
    coverage = _row(
        conn, f'SELECT canonicalization_coverage FROM {SCHEMA}.build_manifest'
    )[0]
    assert isinstance(coverage, dict)
    assert 'by_status' in coverage and 'total' in coverage
    assert coverage['total'] == sum(coverage['by_status'].values())


def test_coverage_counts_are_excluded_from_build_id(conn):
    package_commits, resources, build_id = _row(
        conn,
        f"""
        SELECT package_commits, resources, build_id
        FROM {SCHEMA}.build_manifest
        """,
    )
    # The identity hash covers only the package commits and per-resource content;
    # coverage counts and the translation-table listing must not participate.
    payload = {'package_commits': package_commits, 'resources': resources}
    expected = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(',', ':')).encode()
    ).hexdigest()[:12]
    assert build_id == expected
