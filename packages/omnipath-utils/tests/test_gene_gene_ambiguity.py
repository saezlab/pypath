"""Gene-gene ambiguity in the gene resolver: one identifier, several genes.

The three shapes and their handling (see research R11):
- identical gene copies (same protein, one UniProt, several Entrez genes) — the
  resolver returns EVERY Entrez anchor, so the build can duplicate the mention
  across all of them (no anchor arbitrarily chosen);
- a merged gene — resolves to its surviving gene.

Live checks against a built instance (set ``OMNIPATH_UTILS_TEST_DB``); the
identical-copy examples are stable, well-known human loci.
"""

import os

import pytest

pytestmark = pytest.mark.skipif(
    not os.environ.get('OMNIPATH_UTILS_TEST_DB'),
    reason='set OMNIPATH_UTILS_TEST_DB to a built instance to run',
)


def _session():
    from sqlalchemy.orm import Session

    from omnipath_utils.db._connection import get_engine

    return Session(get_engine(os.environ['OMNIPATH_UTILS_TEST_DB']))


def _anchors(session, accession: str) -> set[int]:
    from sqlalchemy import text

    rows = session.execute(
        text(
            "SELECT DISTINCT entrez FROM omnipath_utils.resolver_gene "
            "WHERE source_type = 'uniprot' AND ncbi_tax_id = 9606 "
            "AND source_id = :ac"
        ),
        {'ac': accession},
    ).all()
    return {int(r[0]) for r in rows}


class TestIdenticalCopies:
    def test_alpha_globin_returns_both_anchors(self):
        # HBA1/HBA2 are identical alpha-globins sharing UniProt P69905; the
        # resolver returns BOTH Entrez genes, deterministically.
        with _session() as session:
            assert _anchors(session, 'P69905') == {3039, 3040}

    def test_hsp70_returns_both_anchors(self):
        # HSPA1A/HSPA1B, identical, share P0DMV8.
        with _session() as session:
            assert _anchors(session, 'P0DMV8') == {3303, 3304}

    def test_copy_loci_are_not_collapsed_to_one(self):
        # The copy shape exists at scale — many human accessions carry exactly
        # two anchors — so the resolver is not silently picking a single gene.
        from sqlalchemy import text

        with _session() as session:
            two_anchor_loci = session.execute(
                text(
                    "SELECT count(*) FROM ("
                    "  SELECT source_id FROM omnipath_utils.resolver_gene"
                    "  WHERE source_type = 'uniprot' AND ncbi_tax_id = 9606"
                    "  GROUP BY source_id HAVING count(DISTINCT entrez) = 2"
                    ") t"
                )
            ).scalar()
        assert two_anchor_loci > 1000


class TestMergedGene:
    def test_merged_gene_maps_to_survivor(self):
        # Entrez 8371 was merged into 19; the history mapping carries the
        # survivor so a merged id can be recovered to the current gene.
        from sqlalchemy import text

        with _session() as session:
            survivor = session.execute(
                text(
                    "SELECT m.target_id FROM omnipath_utils.id_mapping m "
                    "JOIN omnipath_utils.id_type st ON m.source_type_id = st.id "
                    " AND st.name = 'entrez-history' "
                    "JOIN omnipath_utils.id_type tt ON m.target_type_id = tt.id "
                    " AND tt.name = 'entrez' "
                    "WHERE m.ncbi_tax_id = 9606 AND m.source_id = '8371'"
                )
            ).scalar()
        assert survivor == '19'
