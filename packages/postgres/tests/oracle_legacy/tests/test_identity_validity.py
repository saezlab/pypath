"""No canonical identity fails its own value pattern (spec 011 T030, SC-007).

A value like ``InChIKey=none``, ``none``, ``null`` or ``na`` must never become
a canonical identity (spec.md US1 acceptance scenario 4). Two checks prove
it: the registered pattern itself rejects those literal strings (a pure
check, always runnable), and no row in a built database's ``entity`` table
for a chemical actually violates its type's pattern (SC-007, the contract
query in contracts/coverage-acceptance.sql).

Run the DB-backed tests against a built instance::

    DATABASE_URL=postgresql://user:pass@host:5432/dbname \
        uv run --with pytest --with psycopg2-binary \
        pytest tests/test_identity_validity.py -v
"""

from __future__ import annotations

import os
import re

import pytest

DATABASE_URL = os.environ.get('DATABASE_URL')


@pytest.mark.parametrize(
    'garbage_value',
    ['none', 'None', 'NONE', 'null', 'NULL', 'na', 'N/A', 'InChIKey=none'],
)
def test_inchikey_value_pattern_rejects_placeholder_strings(garbage_value):
    from omnipath_utils.mapping._id_types import value_pattern

    pattern = value_pattern('inchikey')
    assert pattern is not None
    assert not re.fullmatch(pattern, garbage_value), (
        f'{garbage_value!r} matches the inchikey value_pattern {pattern!r} -- '
        'a placeholder string must never look like a valid identity'
    )


pg_pytestmark = pytest.mark.skipif(
    not DATABASE_URL,
    reason='DATABASE_URL not set; needs a built database',
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


def _column_exists(conn, table, column):
    return bool(_rows(
        conn,
        """
        SELECT 1 FROM information_schema.columns
        WHERE table_name = %s AND column_name = %s
        """,
        [table, column],
    ))


@pg_pytestmark
def test_vocab_identifier_type_has_value_pattern_column(conn):
    assert _column_exists(conn, 'vocab_identifier_type', 'value_pattern'), (
        'vocab_identifier_type.value_pattern does not exist -- SC-007 (spec '
        '011 T039) needs it to check every canonical identity against its '
        "own namespace's syntax"
    )


@pg_pytestmark
def test_structure_and_chebi_fallback_identity_have_a_registered_value_pattern(conn):
    """SC-007's own query only checks types that *have* a pattern -- a type
    with none is silently skipped, not a failure. WP1 only designates two
    canonical identity types for a chemical: the structure key (InChIKey)
    and the no-structure ChEBI fallback (spec.md US1 acceptance scenario 2).
    Every *other* type a chemical entity's canonical identity might carry
    (bindingdb, foodb, an opaque per-resource accession, ...) is a
    pre-existing fallback this cycle does not touch."""

    from omnipath_build.cv_terms import CHEMICAL_ENTITY_TYPE
    from omnipath_build.duckdb_load import CHEBI_TYPE, STANDARD_INCHI_KEY_TYPE

    if not _column_exists(conn, 'vocab_identifier_type', 'value_pattern'):
        pytest.skip('value_pattern does not exist yet')
    rows = _rows(
        conn,
        """
        SELECT DISTINCT t.name, t.value_pattern
        FROM entity e
        JOIN vocab_identifier_type t ON t.identifier_type_id = e.canonical_identifier_type_id
        JOIN vocab_entity_type et ON et.entity_type_id = e.entity_type_id
        WHERE et.name = %s AND t.name IN (%s, %s)
        """,
        [CHEMICAL_ENTITY_TYPE, CHEBI_TYPE, STANDARD_INCHI_KEY_TYPE],
    )
    assert rows, 'neither structure nor ChEBI-fallback identities found'
    missing = [name for name, pattern in rows if not pattern]
    assert not missing, (
        f'chemical identifier types with no value_pattern: {missing}'
    )


@pg_pytestmark
def test_no_canonical_chemical_identity_violates_its_value_pattern(conn):
    """SC-007's own contract query."""

    from omnipath_build.cv_terms import CHEMICAL_ENTITY_TYPE

    if not _column_exists(conn, 'vocab_identifier_type', 'value_pattern'):
        pytest.skip('value_pattern does not exist yet')
    row = _rows(
        conn,
        """
        SELECT e.canonical_identifier, t.name, t.value_pattern
        FROM entity e
        JOIN vocab_identifier_type t ON t.identifier_type_id = e.canonical_identifier_type_id
        JOIN vocab_entity_type et ON et.entity_type_id = e.entity_type_id
        WHERE et.name = %s
          AND t.value_pattern IS NOT NULL
          AND e.canonical_identifier !~ t.value_pattern
        LIMIT 5
        """,
        [CHEMICAL_ENTITY_TYPE],
    )
    assert not row, f'canonical identities violating their own pattern: {row}'
