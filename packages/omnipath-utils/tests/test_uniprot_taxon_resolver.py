"""The UniProt accession -> organism lookup that lets the network build recover
the species of a protein mention citing an accession but no organism.

Two layers: a contract check that the packaged DDL still exposes the table name
and columns the build reads (a rename would silently break the build's recovery
step), and an optional live-database check (set ``OMNIPATH_UTILS_TEST_DB`` to a
built instance) that known accessions map to their organism.
"""

import importlib.resources as ir
import os

import pytest


def _ddl() -> str:
    return (
        ir.files('omnipath_utils.db.sql')
        .joinpath('resolver_uniprot_taxon.sql')
        .read_text(encoding='utf-8')
    )


class TestContract:
    def test_ddl_is_packaged(self):
        assert 'resolver_uniprot_taxon' in _ddl()

    def test_exposes_accession_and_organism_columns(self):
        ddl = _ddl()
        # The build probes by ``uniprot`` and reads ``ncbi_tax_id`` — renaming
        # either silently disables organism recovery in the network build.
        assert 'AS uniprot' in ddl
        assert 'AS ncbi_tax_id' in ddl

    def test_keyed_by_accession(self):
        assert 'resolver_uniprot_taxon (uniprot)' in _ddl()

    def test_is_a_gene_resolver_ftp_core(self):
        # Built alongside the gene FTP core and gated the same way, so it is
        # (re)built on an FTP reload and reused otherwise.
        from omnipath_utils.db._build import DatabaseBuilder

        source = ir.files('omnipath_utils.db').joinpath('_build.py').read_text(
            encoding='utf-8',
        )
        assert 'sql/resolver_uniprot_taxon.sql' in source
        assert "not present('resolver_uniprot_taxon')" in source
        assert 'gene' in DatabaseBuilder.RESOLVER_NAMES


@pytest.mark.skipif(
    not os.environ.get('OMNIPATH_UTILS_TEST_DB'),
    reason='set OMNIPATH_UTILS_TEST_DB to a built instance to run',
)
class TestLive:
    def _taxon(self, session, accession):
        from sqlalchemy import text

        row = session.execute(
            text(
                'SELECT ncbi_tax_id FROM omnipath_utils.resolver_uniprot_taxon '
                'WHERE uniprot = :ac'
            ),
            {'ac': accession},
        ).fetchone()
        return row[0] if row else None

    def test_known_accessions_carry_their_organism(self):
        from sqlalchemy.orm import Session

        from omnipath_utils.db._connection import get_engine

        engine = get_engine(os.environ['OMNIPATH_UTILS_TEST_DB'])
        with Session(engine) as session:
            # A zebrafish TrEMBL accession and a viral SwissProt accession — both
            # outside the human/mouse core, the exact tail this lookup must cover.
            assert self._taxon(session, 'A0A8M9QHZ6') == 7955
            assert self._taxon(session, 'Q9QBY3') is not None

    def test_one_row_per_accession(self):
        from sqlalchemy import text
        from sqlalchemy.orm import Session

        from omnipath_utils.db._connection import get_engine

        engine = get_engine(os.environ['OMNIPATH_UTILS_TEST_DB'])
        with Session(engine) as session:
            dupes = session.execute(
                text(
                    'SELECT count(*) FROM ('
                    '  SELECT uniprot FROM omnipath_utils.resolver_uniprot_taxon'
                    '  GROUP BY uniprot HAVING count(*) > 1'
                    ') d'
                )
            ).scalar()
        assert dupes == 0


@pytest.mark.skipif(
    not os.environ.get('OMNIPATH_UTILS_TEST_DB'),
    reason='set OMNIPATH_UTILS_TEST_DB to a built instance to run',
)
class TestOrganismQuery:
    def _run(self, accessions):
        from sqlalchemy.orm import Session

        from omnipath_utils.db._connection import get_engine
        from omnipath_utils.db._query import uniprot_organism

        engine = get_engine(os.environ['OMNIPATH_UTILS_TEST_DB'])
        with Session(engine) as session:
            return uniprot_organism(session, accessions)

    def test_known_accessions(self):
        out = self._run(['A0A8M9QHZ6', 'Q9QBY3', 'P04637'])
        assert out['A0A8M9QHZ6'] == 7955
        assert out['Q9QBY3'] == 388906
        assert out['P04637'] == 9606

    def test_isoform_and_curie_are_normalised(self):
        # An isoform suffix and a CURIE prefix both resolve to the base
        # accession's organism, keyed by the input as sent.
        out = self._run(['A0A8M9QHZ6-2', 'uniprot:P04637'])
        assert out['A0A8M9QHZ6-2'] == 7955
        assert out['uniprot:P04637'] == 9606

    def test_unknown_accession_is_none(self):
        out = self._run(['NOTANACCESSION'])
        assert out['NOTANACCESSION'] is None
