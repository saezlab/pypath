"""Tests for the full_uniprot query policy in db._query.translate_ids.

These are unit tests of the curated-first / full-UniProt-fallback orchestration
(FR-025): `_query_table` and `_ftp_table_exists` are monkeypatched so the policy
logic is exercised without a live Postgres.
"""

from collections import defaultdict

import pytest

from omnipath_utils.db import _query


# Curated table resolves TP53; the comprehensive full-UniProt table additionally
# has an extra TP53 accession and a non-human ortholog the curated table lacks.
_CURATED = {'TP53': {'P04637'}}
_FULL = {'TP53': {'P04637', 'X9'}, 'EGFR_CHIMP': {'A1', 'A2'}}


def _fake_query_table(session, table, identifiers, src, tgt, tax):
    data = _CURATED if table.endswith('.id_mapping') else _FULL
    label = 'curated' if data is _CURATED else 'uniprot_ftp'
    res = defaultdict(set)
    backends = set()
    for key, vals in data.items():
        if identifiers is None or key in identifiers:
            res[key] |= vals
            backends.add(label)
    return res, backends


@pytest.fixture(autouse=True)
def _patch(monkeypatch):
    monkeypatch.setattr(_query, '_query_table', _fake_query_table)
    monkeypatch.setattr(_query, '_ftp_table_exists', lambda session: True)
    # genesymbol→uniprot is UniProt-relevant, so the FR-018 gate must allow the
    # FTP fallback these policy tests exercise.
    monkeypatch.setattr(
        _query, '_ftp_types', lambda session: frozenset({'genesymbol', 'uniprot'})
    )


def _translate(ids, **kw):
    return _query.translate_ids(object(), ids, 'genesymbol', 'uniprot', 9606, **kw)


def test_fallback_only_for_unresolved():
    # Default 'fallback': TP53 is resolved by curated, so the full table is NOT
    # consulted for it (it keeps only the curated accession, not the extra X9);
    # the missing ortholog falls back to the full table.
    result, backends = _translate(['TP53', 'EGFR_CHIMP'])
    assert result['TP53'] == {'P04637'}
    assert result['EGFR_CHIMP'] == {'A1', 'A2'}
    assert 'uniprot_ftp' in backends and 'curated' in backends


def test_never_is_curated_only():
    result, backends = _translate(['TP53', 'EGFR_CHIMP'], full_uniprot='never')
    assert result['TP53'] == {'P04637'}
    assert 'EGFR_CHIMP' not in result  # unresolved without the full table
    assert backends == {'curated'}


def test_both_merges_and_dedups():
    # 'both' unions curated + full and deduplicates (set union).
    result, _ = _translate(['TP53'], full_uniprot='both')
    assert result['TP53'] == {'P04637', 'X9'}


def test_only_uses_full_table():
    result, backends = _translate(['TP53'], full_uniprot='only')
    assert result['TP53'] == {'P04637', 'X9'}
    assert backends == {'uniprot_ftp'}


def test_missing_full_table_is_graceful(monkeypatch):
    # When the full-UniProt table is absent, fallback degrades to curated-only.
    monkeypatch.setattr(_query, '_ftp_table_exists', lambda session: False)
    result, _ = _translate(['TP53', 'EGFR_CHIMP'])
    assert result['TP53'] == {'P04637'}
    assert 'EGFR_CHIMP' not in result
