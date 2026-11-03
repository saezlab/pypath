"""T012a -- FR-018: gate the full-UniProt (id_mapping_ftp) fallback.

The ~744 M-row FTP table is 100% UniProt-family, so a pair with a non-FTP type
on either side has zero rows there. The gate skips that fruitless scan and MUST
return results identical to the prior unconditional fallback (R8). The 2,000-id
latency assertion (SC-009) runs only against a configured test DB.
"""

import os
from collections import defaultdict
from unittest.mock import MagicMock

import pytest

import omnipath_utils.db._query as q

# id_types actually present in id_mapping_ftp on the built instance (baseline
# 2026-06-23): only ``uniprot`` as source, mapping to these targets.
FTP_TYPES = frozenset({
    'uniprot', 'chembl', 'drugbank', 'embl', 'embl_id', 'ensg', 'ensp',
    'enst', 'entrez', 'genesymbol', 'gi', 'hgnc', 'kegg', 'refseqn',
    'refseqp', 'uniprot_entry',
})


class TestPredicate:
    def test_name_to_chebi_not_ftp_relevant(self):
        assert q._ftp_relevant('name', 'chebi', FTP_TYPES) is False

    def test_chebi_to_hmdb_not_ftp_relevant(self):
        assert q._ftp_relevant('chebi', 'hmdb', FTP_TYPES) is False

    def test_kegg_to_chebi_not_ftp_relevant(self):
        # kegg IS an FTP type, but chebi is not -> still skip.
        assert q._ftp_relevant('kegg', 'chebi', FTP_TYPES) is False

    def test_genesymbol_to_uniprot_is_ftp_relevant(self):
        assert q._ftp_relevant('genesymbol', 'uniprot', FTP_TYPES) is True

    def test_uniprot_to_entrez_is_ftp_relevant(self):
        assert q._ftp_relevant('uniprot', 'entrez', FTP_TYPES) is True


@pytest.fixture
def record_tables(monkeypatch):
    """Run translate against fake tables; record which tables were queried.

    The fake returns a hit only for the curated ``id_mapping`` so that the
    ``fallback`` path always has 'missing' ids and *would* consult FTP unless
    the gate stops it.
    """
    tables = []

    def fake_query_table(session, table, ids, src, tgt, tax):
        tables.append(table)
        res = defaultdict(set)
        if table.endswith('.id_mapping') and ids:
            # leave at least one id unresolved to trigger the fallback
            for i in ids[1:]:
                res[i].add('HIT')
        return res, {'curated'}

    monkeypatch.setattr(q, '_query_table', fake_query_table)
    monkeypatch.setattr(q, '_ftp_table_exists', lambda s: True)
    monkeypatch.setattr(q, '_ftp_types', lambda s: FTP_TYPES)
    return tables


class TestGateSkipsFtp:
    def test_chebi_to_hmdb_skips_ftp(self, record_tables):
        q.translate_ids(MagicMock(), ['CHEBI:1', 'CHEBI:2'], 'chebi', 'hmdb', 0)
        assert not any('id_mapping_ftp' in t for t in record_tables)

    def test_name_to_chebi_skips_ftp(self, record_tables):
        # routed to id_mapping_long; FTP never involved
        q.translate_ids(MagicMock(), ['Taurine'], 'name', 'chebi', 0)
        assert not any('id_mapping_ftp' in t for t in record_tables)

    def test_genesymbol_to_uniprot_consults_ftp(self, record_tables):
        q.translate_ids(MagicMock(), ['A', 'B'], 'genesymbol', 'uniprot', 9606)
        assert any('id_mapping_ftp' in t for t in record_tables)


class TestResultParity:
    """Gated output equals the unconditional fallback over a matrix (R8).

    FTP holds rows only for (uniprot, X) pairs; for every other pair it returns
    nothing, so gating it off changes no result.
    """

    MATRIX = [
        ('genesymbol', 'uniprot'),   # UniProt-relevant
        ('uniprot', 'entrez'),       # UniProt-relevant
        ('chebi', 'hmdb'),           # chemical <-> chemical
        ('kegg', 'chebi'),           # FTP-type source, non-FTP target
        ('chembl', 'chebi'),         # FTP-type source, non-FTP target
    ]

    def _fake(self, monkeypatch, gate_on):
        def fake_query_table(session, table, ids, src, tgt, tax):
            res = defaultdict(set)
            if ids is None:
                ids = []
            if table.endswith('.id_mapping'):
                for i in ids:
                    res[i].add(f'{src}->{tgt}:curated')
            elif table.endswith('.id_mapping_ftp'):
                # FTP only contains uniprot-family pairs
                if src in FTP_TYPES and tgt in FTP_TYPES:
                    for i in ids:
                        res[i].add(f'{src}->{tgt}:ftp')
            return res, {table.rsplit('.', 1)[-1]}

        monkeypatch.setattr(q, '_query_table', fake_query_table)
        monkeypatch.setattr(q, '_ftp_table_exists', lambda s: True)
        if gate_on:
            monkeypatch.setattr(q, '_ftp_types', lambda s: FTP_TYPES)
        else:
            # disable the gate: every pair is "relevant" -> unconditional
            monkeypatch.setattr(q, '_ftp_relevant', lambda a, b, c: True)
            monkeypatch.setattr(q, '_ftp_types', lambda s: FTP_TYPES)

    @pytest.mark.parametrize('src,tgt', MATRIX)
    def test_gated_equals_unconditional(self, monkeypatch, src, tgt):
        ids = ['X1', 'X2', 'X3']

        self._fake(monkeypatch, gate_on=False)
        ref, _ = q.translate_ids(MagicMock(), list(ids), src, tgt, 0)

        self._fake(monkeypatch, gate_on=True)
        gated, _ = q.translate_ids(MagicMock(), list(ids), src, tgt, 0)

        assert gated == ref


@pytest.mark.skipif(
    not os.environ.get('OMNIPATH_UTILS_TEST_DB'),
    reason='set OMNIPATH_UTILS_TEST_DB to run the SC-009 latency assertion',
)
class TestLatency:
    def test_2000_kegg_to_chebi_under_1s(self):
        import time

        from sqlalchemy.orm import Session

        from omnipath_utils.db._connection import get_engine

        engine = get_engine(os.environ['OMNIPATH_UTILS_TEST_DB'])
        ids = [f'C{i:05d}' for i in range(2000)]
        with Session(engine) as s:
            start = time.time()
            q.translate_ids(s, ids, 'kegg', 'chebi', 0)
            elapsed = time.time() - start
        assert elapsed < 1.0, f'kegg->chebi 2000-batch took {elapsed:.2f}s'
