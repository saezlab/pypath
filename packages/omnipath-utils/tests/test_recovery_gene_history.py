"""T023 -- US2 deprecated-ID recovery: NCBI Gene history (retired Entrez).

A retired Entrez id present in ``gene_history`` recovers to the current id it was
merged into (``recovered=true``, ``recovery_source=gene_history``); a valid id the
primary route already maps never consults recovery (FR-007). DB-free (query layers
monkeypatched).
"""

from collections import defaultdict
from unittest.mock import MagicMock

import pytest

import omnipath_utils.db._query as q


@pytest.fixture
def no_ftp(monkeypatch):
    monkeypatch.setattr(q, '_ftp_table_exists', lambda s: False)
    monkeypatch.setattr(q, '_is_long_query', lambda s, t: False)


def _patch(monkeypatch, primary=None, recovery=None):
    def fake_query_table(session, table, ids, src, tgt, tax):
        res = defaultdict(set)
        for i in (ids or []):
            for v in (primary or {}).get((i, src, tgt), set()):
                res[i].add(v)
        return res, {'curated'} if res else set()

    def fake_recover_query(session, ids, src, tgt, tax):
        res = defaultdict(set)
        for i in ids:
            for v in (recovery or {}).get((i, src, tgt), set()):
                res[i].add(v)
        return res

    monkeypatch.setattr(q, '_query_table', fake_query_table)
    monkeypatch.setattr(q, '_recover_query', fake_recover_query)


def test_retired_entrez_recovers_to_current(monkeypatch, no_ftp):
    # 8371 was discontinued and merged into 19; primary entrez->entrez misses.
    _patch(
        monkeypatch,
        primary={},
        recovery={('8371', 'entrez-history', 'entrez'): {'19'}},
    )
    meta = {}
    result, backends = q.translate_ids(
        MagicMock(), ['8371'], 'entrez', 'entrez', 9606,
        recover=True, recovery_meta=meta,
    )
    assert result['8371'] == {'19'}
    assert meta['8371']['recovered'] is True
    assert meta['8371']['recovery_source'] == 'gene_history'
    assert 'gene_history' in backends


def test_removed_entrez_no_successor_deleted(monkeypatch, no_ftp):
    # A gene removed with no successor is stored self-referentially -> deleted.
    _patch(
        monkeypatch,
        primary={},
        recovery={('999999', 'entrez-history', 'entrez-history'): {'999999'}},
    )
    meta = {}
    result, _ = q.translate_ids(
        MagicMock(), ['999999'], 'entrez', 'entrez', 9606,
        recover=True, recovery_meta=meta,
    )
    assert not result.get('999999')
    assert meta['999999']['deleted'] is True


def test_valid_entrez_does_not_consult_recovery(monkeypatch, no_ftp):
    # Primary route succeeds -> recovery is not consulted, no flags.
    _patch(
        monkeypatch,
        primary={('7157', 'entrez', 'genesymbol'): {'TP53'}},
        recovery={('7157', 'entrez-history', 'entrez'): {'SHOULD_NOT_APPEAR'}},
    )
    meta = {}
    result, _ = q.translate_ids(
        MagicMock(), ['7157'], 'entrez', 'genesymbol', 9606,
        recover=True, recovery_meta=meta,
    )
    assert result['7157'] == {'TP53'}
    assert meta == {}
