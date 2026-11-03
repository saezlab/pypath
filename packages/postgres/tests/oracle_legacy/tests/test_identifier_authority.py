"""The authority relation: who mints which identifier namespace.

Run against a built instance::

    DATABASE_URL=postgresql://omnipath:omnipath@localhost:55432/omnipath \
        uv run --with pytest --with psycopg2-binary \
        pytest tests/test_identifier_authority.py -v
"""

from __future__ import annotations

import os

import pytest

DATABASE_URL = os.environ.get('DATABASE_URL')

pytestmark = pytest.mark.skipif(
    not DATABASE_URL,
    reason='DATABASE_URL not set; the authority relation needs a built database',
)


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
    return _scalar(
        conn,
        "SELECT to_regclass(%s) IS NOT NULL",
        [table],
    )


def test_identifier_authority_table_exists(conn):
    assert _table_exists(conn, 'identifier_authority'), (
        'identifier_authority does not exist. It is populated from the '
        '`mints` declaration on each pypath resource description, collected '
        'in the resource walk (omnipath_build/resources.py).'
    )


def test_identifier_authority_is_populated(conn):
    if not _table_exists(conn, 'identifier_authority'):
        pytest.skip('identifier_authority does not exist yet')
    count = _scalar(conn, 'SELECT count(*) FROM identifier_authority')
    assert count > 0, (
        'identifier_authority exists but is empty. The resource walk has not '
        'populated it -- run the step that collects `mints` declarations.'
    )


def test_one_minting_resource_per_namespace(conn):
    if not _table_exists(conn, 'identifier_authority'):
        pytest.skip('identifier_authority does not exist yet')
    duplicated = _rows(
        conn,
        """
        SELECT identifier_type_id, count(*)
          FROM identifier_authority
         GROUP BY 1
        HAVING count(*) > 1
        """,
    )
    assert duplicated == [], (
        f'{len(duplicated)} namespace(s) have more than one minting '
        f'resource: {duplicated}. A namespace has exactly one authority.'
    )


def test_chebi_mints_chebi(conn):
    if not _table_exists(conn, 'identifier_authority'):
        pytest.skip('identifier_authority does not exist yet')
    authorities = _rows(
        conn,
        """
        SELECT ds.name
          FROM identifier_authority a
          JOIN vocab_identifier_type t USING (identifier_type_id)
          JOIN data_source ds USING (source_id)
         WHERE t.name ILIKE 'chebi%'
        """,
    )
    names = {row[0] for row in authorities}
    assert 'chebi' in names, (
        f'ChEBI does not mint its own namespace. Authorities found: {names}'
    )


def test_kegg_does_not_mint_pubchem_compound(conn):
    """The assertion that would have caught the substance-identifier defect.

    KEGG's `conv/pubchem` endpoint returns PubChem substance identifiers, not
    compound identifiers. KEGG must never be recorded as the minting
    authority for `pubchem_compound`, independent of what the ingest schema
    happens to tag that cross-reference as.
    """
    if not _table_exists(conn, 'identifier_authority'):
        pytest.skip('identifier_authority does not exist yet')
    offenders = _rows(
        conn,
        """
        SELECT ds.name
          FROM identifier_authority a
          JOIN vocab_identifier_type t USING (identifier_type_id)
          JOIN data_source ds USING (source_id)
         WHERE ds.name = 'kegg'
           AND t.name ILIKE 'pubchem compound%'
        """,
    )
    assert offenders == [], (
        'KEGG is recorded as minting pubchem_compound -- it only cites '
        'PubChem substance identifiers via conv/pubchem, and is not that '
        "namespace's authority."
    )
