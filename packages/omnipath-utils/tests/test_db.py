"""Tests for the database module."""

from unittest.mock import MagicMock, patch


class TestDatabaseSchema:
    """Test ORM model definitions."""

    def test_models_importable(self):
        from omnipath_utils.db._schema import (
            IdType,
            IdMapping,
        )

        assert IdType.__tablename__ == 'id_type'
        assert IdMapping.__tablename__ == 'id_mapping'

    def test_idtype_schema(self):
        from omnipath_utils.db._schema import IdType

        cols = {c.name for c in IdType.__table__.columns}
        assert 'name' in cols
        assert 'entity_type' in cols
        assert 'curie_prefix' in cols

    def test_idmapping_indexes(self):
        from omnipath_utils.db._schema import IdMapping

        index_names = {idx.name for idx in IdMapping.__table__.indexes}
        assert 'idx_mapping_lookup' in index_names
        assert 'idx_mapping_reverse' in index_names

    def test_organism_schema(self):
        from omnipath_utils.db._schema import Organism

        cols = {c.name for c in Organism.__table__.columns}
        assert 'ncbi_tax_id' in cols
        assert 'kegg_code' in cols
        assert 'ensembl_name' in cols


class TestConnection:
    def test_get_db_url_default(self):
        from omnipath_utils.db._connection import get_db_url

        url = get_db_url()
        assert 'postgresql' in url

    def test_get_db_url_from_env(self, monkeypatch):
        from omnipath_utils.db._connection import get_db_url

        monkeypatch.setenv(
            'OMNIPATH_UTILS_DB_URL', 'postgresql://test:test@host/db'
        )
        assert get_db_url() == 'postgresql://test:test@host/db'


class TestQuery:
    def test_translate_ids_returns_dict(self):
        """Test translate_ids with a mock session."""
        from omnipath_utils.db._query import translate_ids

        mock_session = MagicMock()
        mock_session.execute.return_value = [
            ('TP53', 'P04637', 'uniprot'),
            ('TP53', 'A0A024R1R8', 'uniprot'),
            ('EGFR', 'P00533', 'uniprot'),
        ]

        result, backends = translate_ids(
            mock_session, ['TP53', 'EGFR'], 'genesymbol', 'uniprot', 9606
        )
        assert 'P04637' in result['TP53']
        assert 'P00533' in result['EGFR']

    def test_translate_ids_empty(self):
        from omnipath_utils.db._query import translate_ids

        # execute() must behave for both the row iteration (empty) and the
        # FTP existence / type probes (.scalar()): the FTP table is absent.
        empty = MagicMock()
        empty.__iter__ = lambda self: iter([])
        empty.scalar.return_value = None
        empty.scalars.return_value.all.return_value = []
        empty.fetchall.return_value = []
        mock_session = MagicMock()
        mock_session.execute.return_value = empty

        result = translate_ids(
            mock_session, ['FAKE'], 'genesymbol', 'uniprot', 9606
        )
        result_dict, backends = result; assert result_dict == {} and backends == set()


class TestDatabaseBuilder:
    @patch('omnipath_utils.db._build.get_engine')
    @patch('omnipath_utils.db._build.ensure_schema')
    def test_builder_creation(self, mock_schema, mock_engine):
        from omnipath_utils.db._build import DatabaseBuilder

        _builder = DatabaseBuilder(db_url='postgresql://test/db')
        mock_engine.assert_called_once()
        mock_schema.assert_called_once()


class TestBuildPipeline:
    @patch('omnipath_utils.db._build.get_engine')
    @patch('omnipath_utils.db._build.ensure_schema')
    def test_populate_from_ftp_importable(self, mock_schema, mock_engine):
        from omnipath_utils.db._build import DatabaseBuilder

        _builder = DatabaseBuilder(db_url='postgresql://test/db')
        assert hasattr(_builder, 'populate_from_ftp')
        assert hasattr(_builder, 'populate_mapping')
        assert hasattr(_builder, 'build_all')
