"""spec 011 T078/C6 -- an InChI (or SMILES) with no stored match still
translates to its structure key, computed live via the chemistry toolkit
rather than requiring it to already be one of our own stored source
values. Degrades (no live computation) without the toolkit -- research R14.
"""

from collections import defaultdict
from unittest.mock import MagicMock

import pytest

import omnipath_utils.db._query as q


@pytest.fixture
def empty_db(monkeypatch):
    """No stored match for anything -- forces the live-compute path."""
    monkeypatch.setattr(
        q, '_query_table',
        lambda session, table, ids, src, tgt, tax: (defaultdict(set), set()),
    )


class TestUnknownInchiStillTranslates:
    def test_unknown_inchi_computes_its_own_key(self, empty_db, monkeypatch):
        from omnipath_utils.mapping import _chemistry

        monkeypatch.setattr(
            _chemistry, 'compute_inchikey',
            lambda **kw: 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N' if kw.get('inchi') else None,
        )
        inchi = 'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3'
        res, _ = q.translate_ids(MagicMock(), [inchi], 'inchi', 'inchikey', 0)
        assert res == {inchi: {'LFQSCWFLJHTTHZ-UHFFFAOYSA-N'}}

    def test_unknown_smiles_computes_its_own_key(self, empty_db, monkeypatch):
        from omnipath_utils.mapping import _chemistry

        monkeypatch.setattr(
            _chemistry, 'compute_inchikey',
            lambda **kw: 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N' if kw.get('smiles') else None,
        )
        res, _ = q.translate_ids(MagicMock(), ['CCO'], 'smiles', 'inchikey', 0)
        assert res == {'CCO': {'LFQSCWFLJHTTHZ-UHFFFAOYSA-N'}}

    def test_only_still_missing_identifiers_get_computed(self, monkeypatch):
        # Taurine's inchi is already stored -- must not be recomputed.
        from omnipath_utils.mapping import _chemistry

        stored = 'InChI=1S/C2H7NO3S/c3-1-2-7(4,5)6/h1-3H2,(H,4,5,6)'
        unknown = 'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3'
        monkeypatch.setattr(
            q, '_query_table',
            lambda session, table, ids, src, tgt, tax: (
                defaultdict(set, {stored: {'XOAAWQZATWQOTB-UHFFFAOYSA-N'}}),
                {'chebi'},
            ),
        )
        calls = []

        def fake_compute(**kw):
            calls.append(kw)
            return 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N'

        monkeypatch.setattr(_chemistry, 'compute_inchikey', fake_compute)
        res, _ = q.translate_ids(
            MagicMock(), [stored, unknown], 'inchi', 'inchikey', 0,
        )
        assert res[stored] == {'XOAAWQZATWQOTB-UHFFFAOYSA-N'}
        assert res[unknown] == {'LFQSCWFLJHTTHZ-UHFFFAOYSA-N'}
        assert len(calls) == 1  # only the unknown one was computed

    def test_garbage_input_computes_nothing(self, empty_db, monkeypatch):
        from omnipath_utils.mapping import _chemistry

        monkeypatch.setattr(_chemistry, 'compute_inchikey', lambda **kw: None)
        res, _ = q.translate_ids(
            MagicMock(), ['not a real inchi'], 'inchi', 'inchikey', 0,
        )
        assert res == {}

    def test_non_inchikey_target_never_triggers_live_compute(
        self, empty_db, monkeypatch,
    ):
        from omnipath_utils.mapping import _chemistry

        calls = []
        monkeypatch.setattr(
            _chemistry, 'compute_inchikey',
            lambda **kw: calls.append(kw) or 'X',
        )
        inchi = 'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3'
        q.translate_ids(MagicMock(), [inchi], 'inchi', 'chebi', 0)
        assert calls == []
