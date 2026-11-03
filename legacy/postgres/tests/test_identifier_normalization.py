"""identifier_evidence carries a normalized form beside the raw one
(spec 011 T028).

WP1's join needs both. The resolver lookup (omnipath-utils, T033) stores its
key already normalized, so the build side must offer the identical
normalized form to join against, while still keeping the verbatim value a
consumer or an audit expects to see. ``value_normalized`` is that second
column. The original ``value`` never changes.

Run against a built instance::

    DATABASE_URL=postgresql://user:pass@host:5432/dbname \
        uv run --with pytest --with psycopg2-binary \
        --with-editable ./omnipath-utils \
        pytest tests/test_identifier_normalization.py -v
"""

from __future__ import annotations

import os

import pytest

DATABASE_URL = os.environ.get('DATABASE_URL')

pytestmark = pytest.mark.skipif(
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


def test_value_normalized_column_exists(conn):
    assert _column_exists(conn, 'identifier_evidence', 'value_normalized'), (
        'identifier_evidence.value_normalized does not exist -- WP1 needs a '
        'normalized-form column beside the raw value (spec 011 T038)'
    )


def test_value_normalized_matches_utils_form_and_value_is_preserved(conn):
    from omnipath_build.duckdb_load import RESOLVER_CHEMICAL_SLUG_TO_IDENTIFIER_TYPE
    from omnipath_utils.mapping._id_types import normalize_identifier

    if not _column_exists(conn, 'identifier_evidence', 'value_normalized'):
        pytest.skip('value_normalized does not exist yet')

    slug_by_type_name = {
        type_name: slug
        for slug, type_name in RESOLVER_CHEMICAL_SLUG_TO_IDENTIFIER_TYPE.items()
    }
    rows = _rows(
        conn,
        """
        SELECT t.name, i.value, i.value_normalized
        FROM identifier_evidence i
        JOIN vocab_identifier_type t ON t.identifier_type_id = i.identifier_type_id
        WHERE t.name = ANY(%s)
        LIMIT 5000
        """,
        [list(slug_by_type_name)],
    )
    assert rows, 'no chemical identifier_evidence rows to check'

    mismatches = []
    for type_name, value, value_normalized in rows:
        expected = normalize_identifier(slug_by_type_name[type_name], value)
        if value_normalized != expected:
            mismatches.append((type_name, value, value_normalized, expected))
        assert value, 'the raw value must never be blanked out'

    assert not mismatches, (
        f'{len(mismatches)} rows have a value_normalized that disagrees with '
        f'omnipath-utils normalize_identifier(), e.g. {mismatches[:5]}'
    )
