"""spec 011 T076 -- contract tests C6-C8 from `contracts/translation-api.md`
(structures as identifier types), run against a built instance (DB-backed,
same convention as test_name_translation.py).
"""

import os

import pytest

pytestmark = pytest.mark.skipif(
    not os.environ.get('OMNIPATH_UTILS_TEST_DB'),
    reason='set OMNIPATH_UTILS_TEST_DB to a built instance to run',
)

rdkit = pytest.importorskip('rdkit', reason='rdkit not installed')


@pytest.fixture(scope='module')
def session():
    from sqlalchemy.orm import Session

    from omnipath_utils.db._connection import get_engine

    engine = get_engine(os.environ['OMNIPATH_UTILS_TEST_DB'])
    with Session(engine) as s:
        yield s


def _translate(session, ids, src, tgt):
    from omnipath_utils.db._query import translate_ids

    res, backends = translate_ids(session, ids, src, tgt, 0)
    return res, backends


def test_c6_an_unknown_inchi_still_translates_to_a_structure_key(session):
    # Ethanol, contracts/translation-api.md's own worked example -- almost
    # certainly not itself a stored source value in id_mapping_long.
    inchi = 'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3'
    res, backends = _translate(session, [inchi], 'inchi', 'inchikey')
    assert res.get(inchi) == {'LFQSCWFLJHTTHZ-UHFFFAOYSA-N'}
    assert 'chemistry_toolkit' in backends


def test_c7_two_equivalent_smiles_spellings_return_the_same_answer(session):
    res_a, _ = _translate(session, ['CCO'], 'smiles', 'inchikey')
    res_b, _ = _translate(session, ['OCC'], 'smiles', 'inchikey')
    assert res_a.get('CCO') == res_b.get('OCC')
    assert res_a.get('CCO')  # both actually answered, not both-empty


def test_c8_a_prefixed_lower_case_structure_key_translates(session):
    # A real chebi<->inchikey pair from the live DB, queried as
    # contracts/translation-api.md's own example shape:
    # 'inchikey=<lowercase key>'.
    from sqlalchemy import text

    row = session.execute(text(
        "SELECT source_id, target_id FROM omnipath_utils.id_mapping m "
        "JOIN omnipath_utils.id_type s ON s.id=m.source_type_id "
        "JOIN omnipath_utils.id_type t ON t.id=m.target_type_id "
        "WHERE s.name='chebi' AND t.name='inchikey' LIMIT 1"
    )).fetchone()
    if not row:
        pytest.skip('no chebi->inchikey rows')
    chebi, inchikey = row
    prefixed = f'inchikey={inchikey.lower()}'
    res, _ = _translate(session, [prefixed], 'inchikey', 'chebi')
    assert chebi in res.get(prefixed, set())


def test_t083_mixed_batch_translates_above_90_percent(session):
    """Verify: a mixed batch of identifiers known to exist in the source
    releases (real ids sampled from id_mapping, one per structure backend)
    translates to inchikey above 90%.
    """
    from sqlalchemy import text

    from omnipath_utils.db._build import DatabaseBuilder

    batch: list[tuple[str, str]] = []  # (id_type, source_id)
    for backend in DatabaseBuilder._STRUCTURE_BACKENDS:
        row = session.execute(text(
            "SELECT s.name, m.source_id FROM omnipath_utils.id_mapping m "
            "JOIN omnipath_utils.id_type s ON s.id=m.source_type_id "
            "JOIN omnipath_utils.id_type t ON t.id=m.target_type_id "
            "WHERE t.name='inchikey' AND s.name=:backend LIMIT 5"
        ), {'backend': backend}).fetchall()
        batch.extend((r[0], r[1]) for r in row)

    if len(batch) < 10:
        pytest.skip(f'too few sampled ids ({len(batch)}) for a meaningful check')

    hits = 0
    for id_type, source_id in batch:
        res, _ = _translate(session, [source_id], id_type, 'inchikey')
        if res.get(source_id):
            hits += 1

    hit_rate = hits / len(batch)
    assert hit_rate >= 0.90, (
        f'{hits}/{len(batch)} = {hit_rate:.1%}, below the 90% target'
    )
