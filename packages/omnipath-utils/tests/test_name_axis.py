"""spec 011 T063/T066 -- the generic 'name' axis unions name, synonym, iupac
and traditional_iupac behind one lookup (FR-016..021).

A name is what the minting resource recommends (the 'preferred' bucket, the
literal 'name' id_type); every other name-like value it records -- including
a systematic IUPAC name -- is a synonym for this purpose (FR-016). These
tests exercise ``_query_name_axis`` directly against a fake ``_query_table``
(no live DB needed, same pattern as ``test_query_routing.py``), and the
response-layer ranking in ``_routes_mapping._rank_name_matches``. The
DB-backed contract tests (C1-C4, C9, C10) live in ``test_name_translation.py``.
"""

from collections import defaultdict
from unittest.mock import MagicMock

import pytest

import omnipath_utils.db._query as q


@pytest.fixture
def fake_db(monkeypatch):
    """A fake ``_query_table`` backed by an explicit (src_type, tgt_type) ->
    {identifier: {targets}} map, so a test can set up exactly which real
    id_type answers which query -- the preferred 'name' bucket, or one of
    the synonym-like types.
    """

    data: dict[tuple[str, str], dict[str, set[str]]] = {}
    calls: list[dict] = []

    def fake_query_table(session, table, ids, src, tgt, tax):
        calls.append({'table': table, 'ids': ids, 'src': src, 'tgt': tgt})
        rows = data.get((src, tgt), {})
        result = defaultdict(set)
        for i in ids or []:
            if i in rows:
                result[i] |= rows[i]
        return result, {f'{src}-backend'} if result else set()

    monkeypatch.setattr(q, '_query_table', fake_query_table)
    return data, calls


class TestForwardNameAxis:
    """Translating FROM a name (FR-017): preferred first, then synonyms."""

    def test_synonym_only_name_answers_through_default_route(self, fake_db):
        # C1: GABA is not a chebi-minted 'name', only a 'synonym'.
        data, _ = fake_db
        data[('synonym', 'chebi')] = {'gaba': {'CHEBI:16865'}}
        result, _, matched = q._query_name_axis(
            MagicMock(), 't', ['gaba'], True, 'chebi', 0, 'default',
        )
        assert result['gaba'] == {'CHEBI:16865'}
        assert matched['gaba'] == 'synonym'

    def test_preferred_scope_does_not_fall_back_to_synonym(self, fake_db):
        # C2: name_scope=preferred must not answer a synonym-only name.
        data, _ = fake_db
        data[('synonym', 'chebi')] = {'gaba': {'CHEBI:16865'}}
        result, _, matched = q._query_name_axis(
            MagicMock(), 't', ['gaba'], True, 'chebi', 0, 'preferred',
        )
        assert result == {}
        assert matched == {}

    def test_preferred_name_does_not_consult_synonym_bucket(self, fake_db):
        # A preferred hit is enough -- no fallback query needed at all.
        data, calls = fake_db
        data[('name', 'chebi')] = {'taurine': {'CHEBI:15891'}}
        result, _, matched = q._query_name_axis(
            MagicMock(), 't', ['taurine'], True, 'chebi', 0, 'default',
        )
        assert result['taurine'] == {'CHEBI:15891'}
        assert matched['taurine'] == 'preferred'
        assert all(c['src'] == 'name' for c in calls)

    def test_matched_as_names_the_answering_bucket(self, fake_db):
        # C3: one preferred hit, one synonym-only hit, in the same call.
        data, _ = fake_db
        data[('name', 'chebi')] = {'taurine': {'CHEBI:15891'}}
        data[('synonym', 'chebi')] = {'gaba': {'CHEBI:16865'}}
        _, _, matched = q._query_name_axis(
            MagicMock(), 't', ['taurine', 'gaba'], True, 'chebi', 0, 'default',
        )
        assert matched == {'taurine': 'preferred', 'gaba': 'synonym'}

    def test_synonym_scope_ignores_preferred_hit(self, fake_db):
        data, _ = fake_db
        data[('name', 'chebi')] = {'taurine': {'CHEBI:15891'}}
        result, _, matched = q._query_name_axis(
            MagicMock(), 't', ['taurine'], True, 'chebi', 0, 'synonym',
        )
        assert result == {}
        assert matched == {}

    def test_any_scope_merges_both_buckets(self, fake_db):
        data, _ = fake_db
        data[('name', 'chebi')] = {'x': {'CHEBI:1'}}
        data[('synonym', 'chebi')] = {'x': {'CHEBI:2'}}
        result, _, matched = q._query_name_axis(
            MagicMock(), 't', ['x'], True, 'chebi', 0, 'any',
        )
        assert result['x'] == {'CHEBI:1', 'CHEBI:2'}
        assert matched['x'] == 'any'

    def test_iupac_and_traditional_iupac_are_synonym_scope(self, fake_db):
        # FR-016: a systematic name is a synonym, not a preferred name.
        data, _ = fake_db
        data[('iupac', 'chebi')] = {'gamma-aminobutanoic acid': {'CHEBI:16865'}}
        result, _, matched = q._query_name_axis(
            MagicMock(), 't', ['gamma-aminobutanoic acid'], True, 'chebi', 0,
            'default',
        )
        assert result['gamma-aminobutanoic acid'] == {'CHEBI:16865'}
        assert matched['gamma-aminobutanoic acid'] == 'synonym'
        result, _, _ = q._query_name_axis(
            MagicMock(), 't', ['gamma-aminobutanoic acid'], True, 'chebi', 0,
            'preferred',
        )
        assert result == {}


class TestReverseNameAxis:
    """Translating TO a name (FR-018): preferred only by default."""

    def test_default_returns_preferred_only_no_synonym_fallback(self, fake_db):
        data, _ = fake_db
        data[('chebi', 'synonym')] = {'CHEBI:16865': {'gaba'}}
        result, _, matched = q._query_name_axis(
            MagicMock(), 't', ['CHEBI:16865'], False, 'chebi', 0, 'default',
        )
        assert result == {}
        assert matched == {}

    def test_synonym_scope_reaches_the_synonym(self, fake_db):
        data, _ = fake_db
        data[('chebi', 'synonym')] = {'CHEBI:16865': {'gaba'}}
        result, _, matched = q._query_name_axis(
            MagicMock(), 't', ['CHEBI:16865'], False, 'chebi', 0, 'synonym',
        )
        assert result['CHEBI:16865'] == {'gaba'}
        assert matched['CHEBI:16865'] == 'synonym'


class TestTranslateIdsNameAxisIntegration:
    """The name axis wired through ``translate_ids`` -- rekeying, tax=0."""

    def test_forward_rekeys_to_original_case(self, fake_db, monkeypatch):
        data, _ = fake_db
        data[('name', 'chebi')] = {'taurine': {'CHEBI:15891'}}
        monkeypatch.setattr(q, '_ftp_types', lambda s: frozenset())
        out, _ = q.translate_ids(MagicMock(), ['Taurine'], 'name', 'chebi', 9606)
        assert out == {'Taurine': {'CHEBI:15891'}}

    def test_name_meta_reports_matched_bucket(self, fake_db, monkeypatch):
        data, _ = fake_db
        data[('synonym', 'chebi')] = {'gaba': {'CHEBI:16865'}}
        monkeypatch.setattr(q, '_ftp_types', lambda s: frozenset())
        meta: dict = {}
        q.translate_ids(
            MagicMock(), ['GABA'], 'name', 'chebi', 0, name_meta=meta,
        )
        assert meta == {'GABA': 'synonym'}


class TestRankNameMatches:
    """Response-layer ranking (FR-021, C9): preferred above synonym."""

    def test_preferred_match_ranks_above_synonym_match(self):
        from omnipath_utils.server._routes_mapping import _rank_name_matches

        mapped = {'gaba': ['CHEBI:16865'], 'taurine': ['CHEBI:15891']}
        matched_as = {'gaba': 'synonym', 'taurine': 'preferred'}
        ranked = _rank_name_matches(mapped, matched_as)
        assert list(ranked) == ['taurine', 'gaba']

    def test_within_a_bucket_fewer_targets_ranks_first(self):
        from omnipath_utils.server._routes_mapping import _rank_name_matches

        mapped = {
            'ambiguous': ['CHEBI:1', 'CHEBI:2'],
            'exact': ['CHEBI:3'],
        }
        matched_as = {'ambiguous': 'synonym', 'exact': 'synonym'}
        ranked = _rank_name_matches(mapped, matched_as)
        assert list(ranked) == ['exact', 'ambiguous']
