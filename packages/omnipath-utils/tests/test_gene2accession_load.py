"""NCBI ``gene2accession`` gene-space coverage — RefSeq RNA/protein to Entrez.

Contract check that the loader and its RefSeq id-types stay wired, plus an
optional live check (``OMNIPATH_UTILS_TEST_DB``) that human RefSeq->Entrez rows
are present and reach the gene anchor.
"""

import os

import pytest


class TestContract:
    def test_loader_is_wired(self):
        from omnipath_utils.db._build import DatabaseBuilder

        assert hasattr(DatabaseBuilder, 'load_gene2accession')

    def test_refseq_id_types_present(self):
        import importlib.resources as ir

        yaml = (
            ir.files('omnipath_utils.data')
            .joinpath('id_types.yaml')
            .read_text(encoding='utf-8')
        )
        assert '\nrefseqn:' in yaml
        assert '\nrefseqp:' in yaml


@pytest.mark.skipif(
    not os.environ.get('OMNIPATH_UTILS_TEST_DB'),
    reason='set OMNIPATH_UTILS_TEST_DB to a built instance to run',
)
class TestLive:
    def _count(self, sql):
        from sqlalchemy import text
        from sqlalchemy.orm import Session

        from omnipath_utils.db._connection import get_engine

        engine = get_engine(os.environ['OMNIPATH_UTILS_TEST_DB'])
        with Session(engine) as session:
            return session.execute(text(sql)).scalar()

    def test_human_refseq_to_entrez_present(self):
        n = self._count(
            'SELECT count(*) FROM omnipath_utils.id_mapping m '
            'JOIN omnipath_utils.id_type st ON m.source_type_id = st.id '
            'JOIN omnipath_utils.id_type tt ON m.target_type_id = tt.id '
            "WHERE st.name IN ('refseqn', 'refseqp') AND tt.name = 'entrez' "
            'AND m.ncbi_tax_id = 9606'
        )
        assert n > 0

    def test_refseq_reaches_gene_anchor(self):
        # A RefSeq accession resolves to an Entrez anchor via resolver_gene.
        n = self._count(
            'SELECT count(*) FROM omnipath_utils.resolver_gene '
            "WHERE source_type IN ('refseqn', 'refseqp') AND ncbi_tax_id = 9606"
        )
        assert n > 0
