"""T041 -- mapping contract (additive): name + structures become mappable.

Exercises the translation core the /mapping/* HTTP layer is built on, against a
built instance (OMNIPATH_UTILS_TEST_DB). Asserts the contract deltas in
contracts/mapping-api.md without needing a running server.
"""

import os

import pytest

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


def test_id_types_table_lists_name_and_structures(session):
    from sqlalchemy import text

    rows = session.execute(text(
        "SELECT name, entity_type FROM omnipath_utils.id_type "
        "WHERE name IN ('name','inchi','smiles')"
    )).fetchall()
    by_name = {r[0]: r[1] for r in rows}
    assert by_name.get('name') == 'small_molecule'
    assert 'inchi' in by_name and 'smiles' in by_name


def test_translate_name_to_chebi(session):
    from omnipath_utils.db._query import translate_ids

    res, _ = translate_ids(session, ['Taurine'], 'name', 'chebi', 0)
    assert res.get('Taurine') == {'CHEBI:15891'}


def test_translate_inchi_to_chebi_roundtrip(session):
    from sqlalchemy import text

    from omnipath_utils.db._query import translate_ids

    row = session.execute(text("""
        SELECT m.source_id, m.target_id FROM omnipath_utils.id_mapping_long m
        JOIN omnipath_utils.id_type s ON s.id=m.source_type_id
        JOIN omnipath_utils.id_type t ON t.id=m.target_type_id
        WHERE s.name='chebi' AND t.name='inchi' LIMIT 1
    """)).fetchone()
    if not row:
        pytest.skip('no chebi->inchi rows')
    chebi, inchi = row
    res, _ = translate_ids(session, [inchi], 'inchi', 'chebi', 0)
    assert chebi in res.get(inchi, set())


def test_identify_name(session):
    from omnipath_utils.db._query import identify_ids

    out = identify_ids(session, ['Taurine'], 0)
    types = {m['id_type'] for m in out.get('Taurine', [])}
    assert 'name' in types or 'synonym' in types
