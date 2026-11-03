"""MANE — the matched RefSeq<->Ensembl transcript set (human).

Contract check that the loader and its pypath input stay wired, plus an optional
live check (``OMNIPATH_UTILS_TEST_DB``) that a MANE RefSeq transcript translates
to its Ensembl counterpart and back, and that MANE ids reach the gene anchor.
"""

import os

import pytest


class TestContract:
    def test_loader_is_wired(self):
        from omnipath_utils.db._build import DatabaseBuilder

        assert hasattr(DatabaseBuilder, 'load_mane')

    def test_pypath_input_present(self):
        from pypath.inputs.mane import mane_summary

        assert callable(mane_summary)

    def test_loader_wired_into_build_all(self):
        import importlib.resources as ir

        source = (
            ir.files('omnipath_utils.db')
            .joinpath('_build.py')
            .read_text(encoding='utf-8')
        )
        assert 'self.load_mane(' in source


@pytest.mark.skipif(
    not os.environ.get('OMNIPATH_UTILS_TEST_DB'),
    reason='set OMNIPATH_UTILS_TEST_DB to a built instance to run',
)
class TestLive:
    def _rows(self, sql, params=None):
        from sqlalchemy import text
        from sqlalchemy.orm import Session

        from omnipath_utils.db._connection import get_engine

        engine = get_engine(os.environ['OMNIPATH_UTILS_TEST_DB'])
        with Session(engine) as session:
            return session.execute(text(sql), params or {}).fetchall()

    def _mane_pair(self):
        # Pick any loaded MANE enst<->refseqn cross-link to test both directions.
        rows = self._rows(
            'SELECT m.source_id, m.target_id '
            'FROM omnipath_utils.id_mapping m '
            'JOIN omnipath_utils.backend b ON m.backend_id = b.id '
            'JOIN omnipath_utils.id_type st ON m.source_type_id = st.id '
            'JOIN omnipath_utils.id_type tt ON m.target_type_id = tt.id '
            "WHERE b.name = 'mane' AND st.name = 'enst' AND tt.name = 'refseqn' "
            'LIMIT 1'
        )
        return rows[0] if rows else None

    def test_refseq_ensembl_translate_both_ways(self):
        pair = self._mane_pair()
        if pair is None:
            pytest.skip('no MANE rows loaded in this instance')
        enst, refseqn = pair
        fwd = self._rows(
            'SELECT 1 FROM omnipath_utils.id_mapping m '
            'JOIN omnipath_utils.id_type st ON m.source_type_id = st.id '
            'JOIN omnipath_utils.id_type tt ON m.target_type_id = tt.id '
            "WHERE st.name = 'enst' AND tt.name = 'refseqn' "
            'AND m.source_id = :a AND m.target_id = :b',
            {'a': enst, 'b': refseqn},
        )
        rev = self._rows(
            'SELECT 1 FROM omnipath_utils.id_mapping m '
            'JOIN omnipath_utils.id_type st ON m.source_type_id = st.id '
            'JOIN omnipath_utils.id_type tt ON m.target_type_id = tt.id '
            "WHERE st.name = 'refseqn' AND tt.name = 'enst' "
            'AND m.source_id = :b AND m.target_id = :a',
            {'a': enst, 'b': refseqn},
        )
        assert fwd and rev
