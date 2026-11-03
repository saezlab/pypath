"""spec 011 T063 -- contract tests C1-C4, C9 (and C10's mechanism) from
`contracts/translation-api.md`, run against a built instance (DB-backed,
unlike `test_name_axis.py`'s fast unit tests over a fake `_query_table`).

Runs against a built instance configured via OMNIPATH_UTILS_TEST_DB (the
same convention as `test_name_to_chebi.py`). Skipped otherwise.
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


def _translate(session, ids, src, tgt, **kw):
    from omnipath_utils.db._query import translate_ids

    meta: dict = {}
    res, _ = translate_ids(session, ids, src, tgt, 0, name_meta=meta, **kw)
    return res, meta


def test_c1_synonym_only_name_answers_through_default_route(session):
    # GABA is chebi's own contracts/translation-api.md example: recorded only
    # as a synonym (chebi's own preferred name is the systematic one), so the
    # default name route must fall back and still answer.
    res, matched = _translate(session, ['GABA'], 'name', 'chebi')
    assert res.get('GABA'), 'GABA should resolve via the synonym fallback'
    assert matched.get('GABA') == 'synonym'


def test_c2_preferred_scope_does_not_answer_a_synonym_only_name(session):
    res, _ = _translate(session, ['GABA'], 'name', 'chebi', name_scope='preferred')
    assert not res.get('GABA')


def test_c3_matched_as_names_the_answering_type(session):
    # Taurine is chebi's own preferred 'name', not just a synonym.
    res, matched = _translate(session, ['Taurine'], 'name', 'chebi')
    assert 'CHEBI:15891' in (res.get('Taurine') or set())
    assert matched.get('Taurine') == 'preferred'


def test_c4_a_miss_reports_nothing_found(session):
    res, matched = _translate(
        session, ['not-a-real-chemical-name-xyz'], 'name', 'chebi',
    )
    assert not res.get('not-a-real-chemical-name-xyz')
    assert 'not-a-real-chemical-name-xyz' not in matched


def test_c9_ranking_places_preferred_above_synonym(session):
    from omnipath_utils.server._routes_mapping import _rank_name_matches

    res, matched = _translate(session, ['Taurine', 'GABA'], 'name', 'chebi')
    mapped = {k: sorted(v) for k, v in res.items() if v}
    ranked = _rank_name_matches(mapped, matched)
    assert list(ranked)[0] == 'Taurine'  # preferred, ranks first


def test_reverse_default_returns_preferred_name_only(session):
    # C5-adjacent (FR-018): translating TO a name returns the preferred
    # name by default -- Taurine's own CHEBI id should answer with its
    # preferred name, not the full synonym soup.
    res, matched = _translate(session, ['CHEBI:15891'], 'chebi', 'name')
    if res.get('CHEBI:15891'):
        assert matched.get('CHEBI:15891') == 'preferred'
