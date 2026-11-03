"""The capability registry: what the build could do, and why not otherwise.

Run against a built instance::

    DATABASE_URL=postgresql://omnipath:omnipath@localhost:55432/omnipath \
        uv run --with pytest --with psycopg2-binary \
        pytest tests/test_build_capability.py -v
"""

from __future__ import annotations

import os

import pytest

DATABASE_URL = os.environ.get('DATABASE_URL')

pytestmark = pytest.mark.skipif(
    not DATABASE_URL,
    reason='DATABASE_URL not set; the capability registry needs a built database',
)

EXPECTED_CAPABILITIES = {
    'lipid_nomenclature',
    'structure_key_computation',
    'structure_substrate',
    'in_database_structure_key',
}


@pytest.fixture(scope='module')
def conn():
    import psycopg2

    connection = psycopg2.connect(DATABASE_URL)
    connection.autocommit = True
    try:
        yield connection
    finally:
        connection.close()


def _rows(conn, query, params=None):
    with conn.cursor() as cur:
        cur.execute(query, params)
        return cur.fetchall()


def _scalar(conn, query, params=None):
    return _rows(conn, query, params)[0][0]


def _table_exists(conn, table):
    return _scalar(conn, 'SELECT to_regclass(%s) IS NOT NULL', [table])


def test_build_capability_table_exists(conn):
    assert _table_exists(conn, 'build_capability'), (
        'build_capability does not exist. It records what the build could '
        'do, so a consumer can tell "no such thing" from "this build could '
        'not compute it" (spec 011 data-model section 8).'
    )


def test_build_capability_records_all_four_capabilities(conn):
    if not _table_exists(conn, 'build_capability'):
        pytest.skip('build_capability does not exist yet')
    rows = _rows(conn, 'SELECT capability FROM build_capability')
    names = {row[0] for row in rows}
    assert names == EXPECTED_CAPABILITIES, (
        f'build_capability has {names}, expected exactly '
        f'{EXPECTED_CAPABILITIES}'
    )


def test_every_row_records_availability_and_provider(conn):
    if not _table_exists(conn, 'build_capability'):
        pytest.skip('build_capability does not exist yet')
    rows = _rows(
        conn,
        'SELECT capability, available, provider, reason FROM build_capability',
    )
    for capability, available, provider, reason in rows:
        assert available in (True, False), (
            f'{capability}: available must be a real boolean, not NULL'
        )
        assert provider, f'{capability}: provider must name who would supply it'
        if not available:
            assert reason, (
                f'{capability} is unavailable but carries no reason -- a '
                'consumer needs to know why, not just that it is absent'
            )


def test_an_unavailable_capability_does_not_fail_the_build(conn):
    """The invariant this table exists to prove: the build completes and
    writes a manifest even when every optional capability is absent.

    Right now every capability genuinely *is* absent on a fresh sandbox
    build (no rdkit-computed structure keys, no lipid nomenclature, no
    metabo phase run) -- the honest state to be in, not a fixture to fake.
    This asserts the build_manifest row exists at all, proving absence
    alone never stopped the build from finishing.
    """
    if not _table_exists(conn, 'build_capability'):
        pytest.skip('build_capability does not exist yet')
    manifest_rows = _scalar(conn, 'SELECT count(*) FROM build_manifest')
    assert manifest_rows == 1, (
        'build_manifest has no row -- a build with absent capabilities '
        'must still complete and write one'
    )


def test_capabilities_mirrored_into_the_manifest(conn):
    if not _table_exists(conn, 'build_capability'):
        pytest.skip('build_capability does not exist yet')
    mirrored = _scalar(conn, 'SELECT capabilities FROM build_manifest')
    assert mirrored is not None, (
        'build_manifest.capabilities is NULL -- the capability registry '
        'must be mirrored into the manifest (data-model section 8)'
    )
    assert len(mirrored) == 4, (
        f'build_manifest.capabilities has {len(mirrored)} entries, expected 4'
    )
