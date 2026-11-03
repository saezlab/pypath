"""T024 -- US2 recovery: Ensembl ID-history, gene-split ambiguity, and the gate.

An old Ensembl gene id recovers to its current id via ``ensembl-history``; a
deprecated id that recovers to **more than one** current identifier returns **all**
candidates flagged ``ambiguous`` (FR-009 clarification — never a silent pick); and
recovery is a gated post-primary fallback (only recoverable source_types, only
still-missing ids). DB-free (query layers monkeypatched).

Note: the ``ensembl-history`` *loader* (pypath bulk stable_id_event input) is a
follow-up; the recovery *query* path already handles ``ensg`` via ``_RECOVERY``, so
these tests exercise it with mocked recovery data.
"""

from collections import defaultdict
from unittest.mock import MagicMock

import pytest

import omnipath_utils.db._query as q


@pytest.fixture
def no_ftp(monkeypatch):
    monkeypatch.setattr(q, '_ftp_table_exists', lambda s: False)
    monkeypatch.setattr(q, '_is_long_query', lambda s, t: False)


def _patch(monkeypatch, primary=None, recovery=None, recorder=None):
    def fake_query_table(session, table, ids, src, tgt, tax):
        res = defaultdict(set)
        for i in (ids or []):
            for v in (primary or {}).get((i, src, tgt), set()):
                res[i].add(v)
        return res, {'curated'} if res else set()

    def fake_recover_query(session, ids, src, tgt, tax):
        if recorder is not None:
            recorder.append((src, tgt))
        res = defaultdict(set)
        for i in ids:
            for v in (recovery or {}).get((i, src, tgt), set()):
                res[i].add(v)
        return res

    monkeypatch.setattr(q, '_query_table', fake_query_table)
    monkeypatch.setattr(q, '_recover_query', fake_recover_query)


def test_old_ensg_recovers_to_current(monkeypatch, no_ftp):
    _patch(
        monkeypatch,
        primary={},
        recovery={('ENSG00000000001', 'ensembl-history', 'ensg'):
                  {'ENSG00000139618'}},
    )
    meta = {}
    result, _ = q.translate_ids(
        MagicMock(), ['ENSG00000000001'], 'ensg', 'ensg', 9606,
        recover=True, recovery_meta=meta,
    )
    assert result['ENSG00000000001'] == {'ENSG00000139618'}
    assert meta['ENSG00000000001']['recovered'] is True
    assert meta['ENSG00000000001']['recovery_source'] == 'ensembl_history'


def test_gene_split_returns_all_candidates_ambiguous(monkeypatch, no_ftp):
    # A retired id that split into two current genes -> both, ambiguous=True.
    _patch(
        monkeypatch,
        primary={},
        recovery={('111', 'entrez-history', 'entrez'): {'222', '333'}},
    )
    meta = {}
    result, _ = q.translate_ids(
        MagicMock(), ['111'], 'entrez', 'entrez', 9606,
        recover=True, recovery_meta=meta,
    )
    assert result['111'] == {'222', '333'}
    assert meta['111']['ambiguous'] is True
    assert meta['111']['recovered'] is True


def test_recovery_gated_to_recoverable_source_types(monkeypatch, no_ftp):
    # 'genesymbol' has no recovery route -> _recover_query never called.
    recorder = []
    _patch(monkeypatch, primary={}, recovery={}, recorder=recorder)
    meta = {}
    result, _ = q.translate_ids(
        MagicMock(), ['NOTATHING'], 'genesymbol', 'entrez', 9606,
        recover=True, recovery_meta=meta,
    )
    assert not result.get('NOTATHING')
    assert meta == {}
    assert recorder == []  # gate: recovery tables not touched for a non-recoverable type


def test_recovery_only_for_missing_ids(monkeypatch, no_ftp):
    # One id maps via the primary route, one is retired: recovery runs only for
    # the missing one.
    recorder = []
    _patch(
        monkeypatch,
        primary={('7157', 'entrez', 'entrez'): {'7157'}},
        recovery={('8371', 'entrez-history', 'entrez'): {'19'}},
        recorder=recorder,
    )
    meta = {}
    result, _ = q.translate_ids(
        MagicMock(), ['7157', '8371'], 'entrez', 'entrez', 9606,
        recover=True, recovery_meta=meta,
    )
    assert result['7157'] == {'7157'}
    assert result['8371'] == {'19'}
    assert set(meta) == {'8371'}  # only the retired id got recovery flags
