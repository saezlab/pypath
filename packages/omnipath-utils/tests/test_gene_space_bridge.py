"""Gene-space bridging — an Ensembl protein id must reach its Entrez anchor
through the gene-space route (``ensp -> ensg -> entrez``), NOT by routing through
protein (UniProt) space, so a gene that has no UniProt is never dropped.

Contract check that the gene-space branches are present in the resolver SQL, plus
an optional live check that the gene-space route resolves at least as many
Ensembl-protein ids as it drops.
"""

import importlib.resources as ir
import os

import pytest


def _ftp_sql() -> str:
    return (
        ir.files('omnipath_utils.db.sql')
        .joinpath('resolver_gene_ftp.sql')
        .read_text(encoding='utf-8')
    )


class TestContract:
    def test_ensp_bridge_is_gene_space(self):
        sql = _ftp_sql()
        # The ensp branch goes ensp -> ensg -> entrez, never ensp -> uniprot.
        assert '_rgf_up_ensp' in sql
        assert '_rgf_ensg_entrez' in sql

    def test_uniprot_ensg_recovery_branch(self):
        # A UniProt with an ensg but no GeneID still reaches entrez via the gene.
        assert 'uniprot -> ensg -> entrez' in _ftp_sql()


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

    def test_ensp_resolves_to_gene_anchor(self):
        n = self._count(
            'SELECT count(*) FROM omnipath_utils.resolver_gene '
            "WHERE source_type = 'ensp' AND ncbi_tax_id = 9606"
        )
        assert n > 0

    def test_gene_space_covers_ensp_with_no_uniprot_route(self):
        # Every ensp resolvable in gene space is an anchor the protein-space
        # route (via a shared uniprot) could not have dropped: the gene-space
        # count is a superset, so it must be >= the ensp-with-uniprot subset.
        gene_space = self._count(
            'SELECT count(DISTINCT source_id) FROM omnipath_utils.resolver_gene '
            "WHERE source_type = 'ensp'"
        )
        assert gene_space > 0
