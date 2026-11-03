"""T022 -- US2 deprecated-ID recovery: UniProt secondary + deleted accessions.

A secondary UniProt AC recovers to its primary (``recovered=true``,
``recovery_source=uniprot_sec``); a deleted AC with no successor is reported as
``deleted`` rather than an unexplained empty result (FR-006/008). The recovery
stage is a post-primary fallback: it never runs for an id the primary route
already mapped, and is off unless ``recover=True`` (FR-007).

DB-free: the primary (``_query_table``) and recovery (``_recover_query``) query
layers are monkeypatched, so the tests exercise the recovery control flow only.
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
    """Install fake primary + recovery query layers."""
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


def test_secondary_uniprot_recovers_to_primary(monkeypatch, no_ftp):
    # Q15086 is a secondary AC of P04637; primary route (uniprot->uniprot) misses.
    _patch(
        monkeypatch,
        primary={},
        recovery={('Q15086', 'uniprot-sec', 'uniprot-pri'): {'P04637'}},
    )
    meta = {}
    result, _ = q.translate_ids(
        MagicMock(), ['Q15086'], 'uniprot', 'uniprot', 9606,
        recover=True, recovery_meta=meta,
    )
    assert result['Q15086'] == {'P04637'}
    assert meta['Q15086']['recovered'] is True
    assert meta['Q15086']['recovery_source'] == 'uniprot_sec'
    assert meta['Q15086']['ambiguous'] is False
    assert meta['Q15086']['deleted'] is False


def test_deleted_uniprot_flagged_deleted(monkeypatch, no_ftp):
    # A deleted AC is stored self-referentially -> reported deleted, not empty.
    _patch(
        monkeypatch,
        primary={},
        recovery={('X99999', 'uniprot-deleted', 'uniprot-deleted'): {'X99999'}},
    )
    meta = {}
    result, _ = q.translate_ids(
        MagicMock(), ['X99999'], 'uniprot', 'uniprot', 9606,
        recover=True, recovery_meta=meta,
    )
    assert not result.get('X99999')
    assert meta['X99999']['deleted'] is True
    assert meta['X99999']['recovered'] is False
    assert meta['X99999']['recovery_source'] == 'uniprot_deleted'


def test_cross_namespace_secondary_to_genesymbol(monkeypatch, no_ftp):
    # secondary -> primary (recovery), then primary -> genesymbol (primary route).
    _patch(
        monkeypatch,
        primary={('P04637', 'uniprot', 'genesymbol'): {'TP53'}},
        recovery={('Q15086', 'uniprot-sec', 'uniprot-pri'): {'P04637'}},
    )
    meta = {}
    result, _ = q.translate_ids(
        MagicMock(), ['Q15086'], 'uniprot', 'genesymbol', 9606,
        recover=True, recovery_meta=meta,
    )
    assert result['Q15086'] == {'TP53'}
    assert meta['Q15086']['recovered'] is True


def test_recovery_off_by_default(monkeypatch, no_ftp):
    _patch(
        monkeypatch,
        primary={},
        recovery={('Q15086', 'uniprot-sec', 'uniprot-pri'): {'P04637'}},
    )
    meta = {}
    result, _ = q.translate_ids(
        MagicMock(), ['Q15086'], 'uniprot', 'uniprot', 9606,
        recovery_meta=meta,  # recover defaults to False
    )
    assert not result.get('Q15086')
    assert meta == {}
