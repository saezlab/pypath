"""NCBI ``gene_info`` gene-space coverage — the authoritative gene symbol (and
synonym) to Entrez map for every organism, independent of UniProt.

A contract check that the loader and its id-type pairs stay wired, and an
optional live check (set ``OMNIPATH_UTILS_TEST_DB`` to a built instance) that
human symbol->Entrez rows are present and reach the gene anchor.
"""

import os

import pytest


class TestContract:
    def test_loader_is_wired(self):
        from omnipath_utils.db._build import DatabaseBuilder

        assert hasattr(DatabaseBuilder, 'load_gene_info')

    def test_id_types_present(self):
        import importlib.resources as ir

        yaml = (
            ir.files('omnipath_utils.data')
            .joinpath('id_types.yaml')
            .read_text(encoding='utf-8')
        )
        assert '\ngenesymbol:' in yaml
        assert '\ngenesymbol-syn:' in yaml
        assert '\nentrez:' in yaml


@pytest.mark.skipif(
    not os.environ.get('OMNIPATH_UTILS_TEST_DB'),
    reason='set OMNIPATH_UTILS_TEST_DB to a built instance to run',
)
class TestLive:
    def _count(self, sql, params=None):
        from sqlalchemy import text
        from sqlalchemy.orm import Session

        from omnipath_utils.db._connection import get_engine

        engine = get_engine(os.environ['OMNIPATH_UTILS_TEST_DB'])
        with Session(engine) as session:
            return session.execute(text(sql), params or {}).scalar()

    def test_human_symbol_to_entrez_present(self):
        n = self._count(
            'SELECT count(*) FROM omnipath_utils.id_mapping m '
            'JOIN omnipath_utils.id_type st ON m.source_type_id = st.id '
            'JOIN omnipath_utils.id_type tt ON m.target_type_id = tt.id '
            "WHERE st.name = 'genesymbol' AND tt.name = 'entrez' "
            'AND m.ncbi_tax_id = 9606'
        )
        assert n > 0

    def test_symbol_reaches_gene_anchor(self):
        # TP53's symbol resolves to its Entrez anchor via resolver_gene.
        n = self._count(
            'SELECT count(*) FROM omnipath_utils.resolver_gene '
            "WHERE source_type = 'genesymbol' AND source_id = 'TP53' "
            "AND ncbi_tax_id = 9606 AND entrez = '7157'"
        )
        assert n > 0
