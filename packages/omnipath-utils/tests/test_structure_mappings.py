"""C9/T027b -- structures (InChI/SMILES) loaded & stored verbatim (DB-backed).

Runs against a built instance via OMNIPATH_UTILS_TEST_DB.
"""

import os

import pytest
from sqlalchemy import text

pytestmark = pytest.mark.skipif(
    not os.environ.get('OMNIPATH_UTILS_TEST_DB'),
    reason='set OMNIPATH_UTILS_TEST_DB to a built instance to run',
)


@pytest.fixture(scope='module')
def session():
    from sqlalchemy.orm import Session

    from omnipath_utils.db._connection import get_engine

    engine = get_engine(os.environ['OMNIPATH_UTILS_TEST_DB'])
    with Session(engine) as s:
        yield s


def test_inchi_and_smiles_have_rows(session):
    rows = session.execute(text("""
        SELECT s.name AS src, t.name AS tgt, count(*) n
        FROM omnipath_utils.id_mapping_long m
        JOIN omnipath_utils.id_type s ON s.id = m.source_type_id
        JOIN omnipath_utils.id_type t ON t.id = m.target_type_id
        WHERE 'inchi' IN (s.name, t.name) OR 'smiles' IN (s.name, t.name)
        GROUP BY 1, 2
    """)).fetchall()
    pairs = {(r[0], r[1]): r[2] for r in rows}
    # inchi<->chebi and smiles<->chebi must exist (at least one direction each)
    assert pairs.get(('chebi', 'inchi'), 0) > 0 or pairs.get(('inchi', 'chebi'), 0) > 0
    assert pairs.get(('chebi', 'smiles'), 0) > 0 or pairs.get(('smiles', 'chebi'), 0) > 0


def test_structures_stored_verbatim(session):
    """A lower()/upper() round-trip would change a structure -> case preserved."""
    rows = session.execute(text("""
        SELECT source_id FROM omnipath_utils.id_mapping_long m
        JOIN omnipath_utils.id_type s ON s.id = m.source_type_id
        WHERE s.name IN ('inchi','smiles') LIMIT 500
    """)).fetchall()
    assert rows, 'no structure source rows'
    # at least one structure must contain mixed case (proving no folding)
    assert any(
        v[0] != v[0].lower() and v[0] != v[0].upper() for v in rows
    ), 'structure keys appear case-folded'


def test_inchi_resolves_to_chebi(session):
    from omnipath_utils.db._query import translate_ids

    # take a real chebi->inchi pair, then resolve the inchi back to chebi
    row = session.execute(text("""
        SELECT m.source_id, m.target_id FROM omnipath_utils.id_mapping_long m
        JOIN omnipath_utils.id_type s ON s.id = m.source_type_id
        JOIN omnipath_utils.id_type t ON t.id = m.target_type_id
        WHERE s.name='chebi' AND t.name='inchi' LIMIT 1
    """)).fetchone()
    if not row:
        pytest.skip('no chebi->inchi rows')
    chebi, inchi = row
    res, _ = translate_ids(session, [inchi], 'inchi', 'chebi', 0)
    assert chebi in res.get(inchi, set())
