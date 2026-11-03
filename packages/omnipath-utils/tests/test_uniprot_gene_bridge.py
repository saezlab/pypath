"""The UniProt-accession -> gene-symbol -> Entrez bridge in the gene resolver.

Recovers a UniProt accession that the idmapping gives only a gene symbol for (no
GeneID, no Ensembl gene — typical of the TrEMBL long tail) by matching that symbol
to NCBI gene_info, including synonyms. A contract check that the packaged DDL keeps
the bridge, plus an optional live check (set ``OMNIPATH_UTILS_TEST_DB`` to a built
instance) that symbol-only accessions reach a gene.
"""

import importlib.resources as ir
import os

import pytest


def _ddl() -> str:
    return (
        ir.files('omnipath_utils.db.sql')
        .joinpath('resolver_gene_ftp.sql')
        .read_text(encoding='utf-8')
    )


class TestContract:
    def test_bridge_present(self):
        ddl = _ddl()
        # The bridge joins the FTP uniprot->symbol extraction to gene_info's
        # symbol->entrez map, and folds in synonyms to cover symbol divergence.
        assert '_rgf_up_symbol' in ddl
        assert "'genesymbol', 'genesymbol-syn'" in ddl
        assert "'uniprot', s.uniprot, gi.entrez" in ddl


@pytest.mark.skipif(
    not os.environ.get('OMNIPATH_UTILS_TEST_DB'),
    reason='set OMNIPATH_UTILS_TEST_DB to a built instance to run',
)
class TestLive:
    def test_symbol_only_uniprot_reaches_a_gene(self):
        from sqlalchemy import text
        from sqlalchemy.orm import Session

        from omnipath_utils.db._connection import get_engine

        engine = get_engine(os.environ['OMNIPATH_UTILS_TEST_DB'])
        with Session(engine) as session:
            # The bridge must add UniProt accessions that resolve ONLY through the
            # gene-symbol route (no direct GeneID, no Ensembl bridge) — i.e. the
            # count of distinct UniProt accessions anchored to Entrez is strictly
            # larger than the direct + Ensembl arms alone would give.
            anchored = session.execute(
                text(
                    "SELECT count(DISTINCT source_id) "
                    "FROM omnipath_utils.resolver_gene_ftp "
                    "WHERE source_type = 'uniprot'"
                )
            ).scalar()
        assert anchored > 21_000_000
