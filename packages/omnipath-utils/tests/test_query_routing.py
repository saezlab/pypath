"""T011 -- query routing between id_mapping and id_mapping_long.

Name/structure queries read the long-value table; database-ID queries keep
reading id_mapping. The query term is lowercased only for name source types.
"""

from collections import defaultdict
from unittest.mock import MagicMock

import pytest

import omnipath_utils.db._query as q


@pytest.fixture
def spy(monkeypatch):
    """Record every ``_query_table`` call as ``(table, ids, tax)``."""
    calls = []

    def fake_query_table(session, table, ids, src, tgt, tax):
        calls.append({'table': table, 'ids': ids, 'src': src,
                      'tgt': tgt, 'tax': tax})
        return defaultdict(set), set()

    monkeypatch.setattr(q, '_query_table', fake_query_table)
    # Keep the FTP gate off the DB in unit tests.
    monkeypatch.setattr(q, '_ftp_types', lambda s: frozenset())
    return calls


class TestLongRouting:
    def test_name_to_chebi_reads_long_lowercased(self, spy):
        # spec 011 FR-017: 'name' is the generic axis -- preferred first (one
        # call), falling back to the synonym/iupac/traditional_iupac bucket
        # (three more) since the spy answers nothing.
        q.translate_ids(MagicMock(), ['Taurine'], 'name', 'chebi', 9606)
        assert len(spy) == 4
        assert all(c['table'].endswith('id_mapping_long') for c in spy)
        assert all(c['ids'] == ['taurine'] for c in spy)  # lowercased (FR-002)
        assert all(c['tax'] == 0 for c in spy)  # organism-agnostic
        assert spy[0]['src'] == 'name'  # the preferred bucket is tried first

    def test_synonym_to_chebi_reads_long(self, spy):
        q.translate_ids(MagicMock(), ['ATP'], 'synonym', 'chebi', 9606)
        assert spy[0]['table'].endswith('id_mapping_long')
        assert spy[0]['ids'] == ['atp']

    def test_inchi_to_chebi_reads_long_verbatim(self, spy):
        inchi = 'InChI=1S/H2O/h1H2'
        q.translate_ids(MagicMock(), [inchi], 'inchi', 'chebi', 0)
        assert spy[0]['table'].endswith('id_mapping_long')
        assert spy[0]['ids'] == [inchi]          # NOT lowercased (FR-019)

    def test_chebi_to_inchi_reads_long(self, spy):
        q.translate_ids(MagicMock(), ['CHEBI:15377'], 'chebi', 'inchi', 0)
        assert spy[0]['table'].endswith('id_mapping_long')
        # chebi is not a name type -> verbatim
        assert spy[0]['ids'] == ['CHEBI:15377']

    def test_chebi_to_name_not_lowercased(self, spy):
        q.translate_ids(MagicMock(), ['CHEBI:15377'], 'chebi', 'name', 0)
        assert spy[0]['table'].endswith('id_mapping_long')
        assert spy[0]['ids'] == ['CHEBI:15377']


class TestNameResultRekeying:
    def test_response_keyed_by_original_input(self, monkeypatch):
        # the long table stores/returns the lowercased key; the response must
        # be re-keyed back to the caller's original-case identifier.
        def fake_query_table(session, table, ids, src, tgt, tax):
            res = defaultdict(set)
            res['taurine'].add('CHEBI:15891')
            return res, {'chebi'}

        monkeypatch.setattr(q, '_query_table', fake_query_table)
        for variant in ('Taurine', 'taurine', 'TAURINE', ' Taurine '):
            out, _ = q.translate_ids(MagicMock(), [variant], 'name', 'chebi', 0)
            assert out == {variant: {'CHEBI:15891'}}

    def test_structure_keys_not_rekeyed(self, monkeypatch):
        inchi = 'InChI=1S/H2O/h1H2'

        def fake_query_table(session, table, ids, src, tgt, tax):
            res = defaultdict(set)
            res[inchi].add('CHEBI:15377')
            return res, {'chebi'}

        monkeypatch.setattr(q, '_query_table', fake_query_table)
        out, _ = q.translate_ids(MagicMock(), [inchi], 'inchi', 'chebi', 0)
        assert out == {inchi: {'CHEBI:15377'}}


class TestDatabaseIdRouting:
    def test_chebi_to_hmdb_reads_id_mapping(self, spy):
        q.translate_ids(
            MagicMock(), ['CHEBI:15377'], 'chebi', 'hmdb', 0,
            full_uniprot='never',
        )
        assert len(spy) == 1
        assert spy[0]['table'].endswith('.id_mapping')
        assert not spy[0]['table'].endswith('id_mapping_long')

    def test_database_id_query_never_touches_long(self, spy):
        q.translate_ids(
            MagicMock(), ['P00533'], 'genesymbol', 'uniprot', 9606,
            full_uniprot='never',
        )
        tables = [c['table'] for c in spy]
        assert all('id_mapping_long' not in t for t in tables)
