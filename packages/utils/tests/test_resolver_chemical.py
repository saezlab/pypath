"""The chemical resolver projection normalizes and deduplicates its keys
(spec 011 T027).

WP1's whole join depends on this. omnipath-build looks up a source
identifier in ``resolver_chemical`` using the *normalized* form (spec 011
T024-T026), so a row stored under the raw, un-normalized form never
matches. And a triple repeated under UNION ALL wastes the join without
changing its result -- the projection should carry each (source_type,
source_id, inchikey) once.

Runs against a built instance via OMNIPATH_UTILS_TEST_DB.
"""

from __future__ import annotations

import os

import pytest
from sqlalchemy import text

from omnipath_utils.mapping._id_types import value_pattern

pytestmark = pytest.mark.skipif(
    not os.environ.get('OMNIPATH_UTILS_TEST_DB'),
    reason='set OMNIPATH_UTILS_TEST_DB to a built instance to run',
)

# Namespaces with a registered value_pattern (spec 011 T024-T026) that
# resolver_chemical actually carries source_ids for.
_CHECKED_TYPES = ['chebi', 'hmdb', 'lipidmaps', 'kegg', 'pubchem', 'chembl']


@pytest.fixture(scope='module')
def session():
    from sqlalchemy.orm import Session

    from omnipath_utils.db._connection import get_engine

    engine = get_engine(os.environ['OMNIPATH_UTILS_TEST_DB'])
    with Session(engine) as s:
        yield s


@pytest.mark.parametrize('source_type', _CHECKED_TYPES)
def test_source_id_is_normalized(session, source_type):
    """Every stored source_id already matches its own value_pattern.

    Checked in SQL, not by pulling every row through Python's
    normalize_identifier -- resolver_chemical carries over 120M rows for
    pubchem alone.
    """

    pattern = value_pattern(source_type)
    assert pattern is not None, f'{source_type} has no registered value_pattern'

    row = session.execute(
        text("""
            SELECT source_id FROM omnipath_utils.resolver_chemical
            WHERE source_type = :source_type
              AND source_id !~ :pattern
            LIMIT 5
        """),
        {'source_type': source_type, 'pattern': f'^(?:{pattern})$'},
    ).fetchall()
    assert not row, (
        f'{source_type}: resolver_chemical source_id values not matching '
        f'value_pattern {pattern!r}, e.g. {[r[0] for r in row]}'
    )


def test_no_duplicate_triples(session):
    """UNION ALL across direct/bridge_pubchem/bridge_chebi must not repeat a
    (source_type, source_id, inchikey) triple."""

    row = session.execute(text("""
        SELECT source_type, source_id, inchikey, count(*) AS n
        FROM omnipath_utils.resolver_chemical
        GROUP BY 1, 2, 3
        HAVING count(*) > 1
        ORDER BY n DESC
        LIMIT 1
    """)).fetchone()
    assert row is None, f'duplicate resolver_chemical triple: {tuple(row)}'
