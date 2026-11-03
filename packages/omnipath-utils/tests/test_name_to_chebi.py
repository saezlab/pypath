"""T013 / C4 -- US1 name->ChEBI integration (DB-backed).

Runs against a built instance configured via OMNIPATH_UTILS_TEST_DB (the utils2
build). Skipped otherwise. Asserts case-insensitive name->chebi resolution and
that a many->one synonym returns the complete candidate set (FR-014).
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


def _translate(session, ids, src, tgt):
    from omnipath_utils.db._query import translate_ids

    res, _ = translate_ids(session, ids, src, tgt, 0)
    return res


def test_taurine_name_to_chebi(session):
    res = _translate(session, ['Taurine'], 'name', 'chebi')
    assert res.get('Taurine') == {'CHEBI:15891'}


def test_case_insensitive(session):
    # C4/SC-007: lower/upper/mixed map identically, keyed by the original input
    out = {}
    for v in ('Taurine', 'taurine', 'TAURINE'):
        out[v] = _translate(session, [v], 'name', 'chebi').get(v)
    assert out['Taurine'] == out['taurine'] == out['TAURINE'] == {'CHEBI:15891'}


def test_synonym_to_chebi(session):
    res = _translate(session, ['ATP'], 'synonym', 'chebi')
    assert res.get('ATP'), 'ATP synonym should resolve to >=1 ChEBI'


def test_many_to_one_synonym_full_candidate_set(session):
    # FR-014: a synonym shared by >=2 ChEBI returns ALL candidates, not one.
    res = _translate(session, ['ATP'], 'synonym', 'chebi')
    assert len(res.get('ATP', set())) >= 2


def test_sample_resolution_at_least_80pct(session):
    """C6/SC-001: >=80% of the 324-name sample resolves to >=1 ChEBI."""
    import json

    sample_path = os.environ.get('OMNIPATH_UTILS_SAMPLE_NAMES')
    if not sample_path or not os.path.exists(sample_path):
        pytest.skip('set OMNIPATH_UTILS_SAMPLE_NAMES to the 324-name JSON')
    names = json.load(open(sample_path))
    resolved = set()
    for src in ('name', 'synonym', 'iupac', 'traditional_iupac'):
        res = _translate(session, names, src, 'chebi')
        resolved.update(n for n, hits in res.items() if hits)
    assert len(resolved) / len(names) >= 0.80
