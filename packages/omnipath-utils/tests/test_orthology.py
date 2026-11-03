"""Tests for the orthology module."""

from unittest.mock import MagicMock, patch

from omnipath_utils.orthology._manager import OrthologyTable, OrthologyManager


class TestOrthologyTable:
    def test_create(self):
        table = OrthologyTable(
            data={'TP53': {'Trp53'}},
            source=9606,
            target=10090,
            id_type='genesymbol',
            resource='hcop',
        )
        assert len(table) == 1
        assert table['TP53'] == {'Trp53'}
        assert table['FAKE'] == set()


class TestOrthologyManager:
    def setup_method(self):
        OrthologyManager._instance = None

    def test_singleton(self):
        m1 = OrthologyManager.get()
        m2 = OrthologyManager.get()
        assert m1 is m2

    def test_translate_with_preloaded_table(self):
        mgr = OrthologyManager()
        table = OrthologyTable(
            data={'TP53': {'Trp53'}, 'EGFR': {'Egfr'}},
            source=9606,
            target=10090,
            id_type='genesymbol',
            resource='hcop',
        )
        mgr._tables[(9606, 10090, 'genesymbol', 'hcop')] = table

        result = mgr.translate(
            ['TP53', 'EGFR', 'FAKE'], source=9606, target=10090
        )
        assert result['TP53'] == {'Trp53'}
        assert result['EGFR'] == {'Egfr'}
        assert result['FAKE'] == set()

    def test_translate_with_metadata_filter(self):
        mgr = OrthologyManager()
        meta = {
            'TP53': {
                'Trp53': {'rel_type': '1:1', 'n_sources': 10},
                'Trp53-like': {'rel_type': '1:n', 'n_sources': 2},
            },
        }
        table = OrthologyTable(
            data={'TP53': {'Trp53', 'Trp53-like'}},
            source=9606,
            target=10090,
            id_type='genesymbol',
            resource='hcop',
            metadata=meta,
        )
        mgr._tables[(9606, 10090, 'genesymbol', 'hcop')] = table

        # Filter to 1:1 only
        result = mgr.translate(
            ['TP53'],
            source=9606,
            target=10090,
            rel_type={'1:1'},
        )
        assert result['TP53'] == {'Trp53'}

    @patch.object(OrthologyManager, '_load_table')
    def test_resource_selection(self, mock_load):
        mock_load.return_value = None
        mgr = OrthologyManager()

        mgr.translate(['TP53'], source=9606, target=10090, resource='oma')

        # Should only try 'oma', not all resources
        calls = mock_load.call_args_list
        assert len(calls) == 1
        assert (
            calls[0].args[3] == 'oma'
            or calls[0].kwargs.get('resource') == 'oma'
        )


class TestOrthologyPublicAPI:
    def setup_method(self):
        OrthologyManager._instance = None

    @patch.object(OrthologyManager, '_get_table')
    def test_translate_function(self, mock_table):
        table = OrthologyTable(
            data={'TP53': {'Trp53'}},
            source=9606,
            target=10090,
            id_type='genesymbol',
            resource='hcop',
        )
        mock_table.return_value = table

        from omnipath_utils.orthology import translate

        result = translate(['TP53'], source=9606, target=10090)
        assert result['TP53'] == {'Trp53'}


class TestOrthologySchema:
    def test_orthology_model_importable(self):
        from omnipath_utils.db._schema import Orthology

        assert Orthology.__tablename__ == 'orthology'
        cols = {c.name for c in Orthology.__table__.columns}
        assert 'source_id' in cols
        assert 'target_id' in cols
        assert 'source_tax_id' in cols
        assert 'rel_type' in cols
        assert 'n_sources' in cols
        assert 'support' in cols


class TestOrthologyBackends:
    def setup_method(self):
        OrthologyManager._instance = None

    @patch.object(OrthologyManager, '_load_table')
    def test_default_resource_order(self, mock_load):
        mock_load.return_value = None
        mgr = OrthologyManager()
        mgr.translate(['TP53'], source=9606, target=10090)

        # Should try all resources in order
        resources_tried = [c.args[3] for c in mock_load.call_args_list]
        assert resources_tried == [
            'hcop',
            'ensembl',
            'oma',
            'orthodb',
            'alliance',
            'homologene',
        ]

    def test_orthodb_backend_exists(self):
        mgr = OrthologyManager()
        assert hasattr(mgr, '_load_orthodb')
        assert callable(mgr._load_orthodb)

    def test_alliance_backend_exists(self):
        mgr = OrthologyManager()
        assert hasattr(mgr, '_load_alliance')
        assert callable(mgr._load_alliance)

    def test_orthodb_loads_table(self):
        mock_mod = MagicMock()
        mock_mod.orthodb_orthologs.return_value = {'TP53': {'Trp53'}}

        with patch.dict(
            'sys.modules',
            {
                'pypath': MagicMock(),
                'pypath.inputs': MagicMock(),
                'pypath.inputs.orthodb': mock_mod,
            },
        ):
            mgr = OrthologyManager()
            table = mgr._load_orthodb(9606, 10090, 'genesymbol')
            assert table is not None
            assert table.resource == 'orthodb'
            assert table['TP53'] == {'Trp53'}
            mock_mod.orthodb_orthologs.assert_called_once_with(
                source=9606,
                target=10090,
                id_type='genesymbol',
            )

    def test_orthodb_empty_returns_none(self):
        mock_mod = MagicMock()
        mock_mod.orthodb_orthologs.return_value = {}

        with patch.dict(
            'sys.modules',
            {
                'pypath': MagicMock(),
                'pypath.inputs': MagicMock(),
                'pypath.inputs.orthodb': mock_mod,
            },
        ):
            mgr = OrthologyManager()
            table = mgr._load_orthodb(9606, 10090, 'genesymbol')
            assert table is None

    def test_alliance_loads_table(self):
        mock_mod = MagicMock()
        mock_mod.alliance_dict.return_value = {'TP53': {'Trp53'}}

        with patch.dict(
            'sys.modules',
            {
                'pypath': MagicMock(),
                'pypath.inputs': MagicMock(),
                'pypath.inputs.alliance': mock_mod,
            },
        ):
            mgr = OrthologyManager()
            table = mgr._load_alliance(9606, 10090, 'genesymbol')
            assert table is not None
            assert table.resource == 'alliance'
            assert table['TP53'] == {'Trp53'}
            mock_mod.alliance_dict.assert_called_once_with(
                source=9606,
                target=10090,
                id_type='genesymbol',
            )

    def test_alliance_empty_returns_none(self):
        mock_mod = MagicMock()
        mock_mod.alliance_dict.return_value = {}

        with patch.dict(
            'sys.modules',
            {
                'pypath': MagicMock(),
                'pypath.inputs': MagicMock(),
                'pypath.inputs.alliance': mock_mod,
            },
        ):
            mgr = OrthologyManager()
            table = mgr._load_alliance(9606, 10090, 'genesymbol')
            assert table is None
