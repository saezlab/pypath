"""Tests for the unified translate_core function and new API parameters."""

from unittest.mock import MagicMock, patch

from omnipath_utils.mapping._table import MappingTable
from omnipath_utils.mapping._mapper import Mapper
from omnipath_utils.mapping._translate import translate_core


class TestTranslateCore:
    """Tests for the unified translate_core function."""

    def setup_method(self):
        """Reset singleton for test isolation."""
        Mapper._instance = None

    @patch('omnipath_utils.mapping._mapper.Mapper._load_table')
    def test_basic_batch_translation(self, mock_load):
        mock_load.return_value = None
        mapper = Mapper()

        table = MappingTable(
            data={'P04637': {'TP53'}, 'P00533': {'EGFR'}},
            id_type='uniprot',
            target_id_type='genesymbol',
            ncbi_tax_id=9606,
        )
        mapper.tables[table.key] = table
        Mapper._instance = mapper

        result = translate_core(
            ['P04637', 'P00533'],
            'uniprot',
            'genesymbol',
        )
        assert result['P04637'] == {'TP53'}
        assert result['P00533'] == {'EGFR'}

    @patch('omnipath_utils.mapping._mapper.Mapper._load_table')
    def test_raw_mode_skips_fallbacks(self, mock_load):
        mock_load.return_value = None
        mapper = Mapper()

        table = MappingTable(
            data={'TP53': {'P04637'}},
            id_type='genesymbol',
            target_id_type='uniprot',
            ncbi_tax_id=9606,
        )
        mapper.tables[table.key] = table
        Mapper._instance = mapper

        # Raw mode: "tp53" (lowercase) should NOT trigger case fallback
        result = translate_core(
            ['tp53'],
            'genesymbol',
            'uniprot',
            raw=True,
        )
        # Direct lookup for "tp53" in the table returns empty
        assert result['tp53'] == set()

    @patch('omnipath_utils.mapping._mapper.Mapper._load_table')
    def test_raw_mode_direct_hit(self, mock_load):
        mock_load.return_value = None
        mapper = Mapper()

        table = MappingTable(
            data={'TP53': {'P04637'}},
            id_type='genesymbol',
            target_id_type='uniprot',
            ncbi_tax_id=9606,
        )
        mapper.tables[table.key] = table
        Mapper._instance = mapper

        # Raw mode: exact match should still work
        result = translate_core(
            ['TP53'],
            'genesymbol',
            'uniprot',
            raw=True,
        )
        assert result['TP53'] == {'P04637'}

    @patch('omnipath_utils.mapping._mapper.Mapper._load_table')
    def test_raw_no_table_returns_empty(self, mock_load):
        mock_load.return_value = None
        mapper = Mapper()
        Mapper._instance = mapper

        result = translate_core(
            ['TP53'],
            'genesymbol',
            'uniprot',
            raw=True,
        )
        assert result['TP53'] == set()

    @patch('omnipath_utils.mapping._mapper.Mapper._load_table')
    def test_fallback_for_missing_ids(self, mock_load):
        mock_load.return_value = None
        mapper = Mapper()

        table = MappingTable(
            data={'TP53': {'P04637'}},
            id_type='genesymbol',
            target_id_type='uniprot',
            ncbi_tax_id=9606,
        )
        mapper.tables[table.key] = table
        Mapper._instance = mapper

        # "EGFR" is not in the table, so it should go through
        # map_name fallback. Since no other tables are loaded,
        # it will return empty.
        result = translate_core(
            ['TP53', 'EGFR'],
            'genesymbol',
            'uniprot',
        )
        assert result['TP53'] == {'P04637'}
        # EGFR goes through map_name which won't find it either
        assert 'EGFR' in result

    @patch('omnipath_utils.mapping._mapper.Mapper._load_table')
    def test_vectorized_hits_skip_fallback(self, mock_load):
        mock_load.return_value = None
        mapper = Mapper()

        table = MappingTable(
            data={
                'P04637': {'TP53'},
                'P00533': {'EGFR'},
                'P38398': {'BRCA1'},
            },
            id_type='uniprot',
            target_id_type='genesymbol',
            ncbi_tax_id=9606,
        )
        mapper.tables[table.key] = table
        Mapper._instance = mapper

        # All IDs are in the table, no fallback needed
        result = translate_core(
            ['P04637', 'P00533', 'P38398'],
            'uniprot',
            'genesymbol',
        )
        assert result['P04637'] == {'TP53'}
        assert result['P00533'] == {'EGFR'}
        assert result['P38398'] == {'BRCA1'}


class TestWhichTableBackend:
    """Test the backend parameter on which_table."""

    def setup_method(self):
        Mapper._instance = None

    @patch('omnipath_utils.mapping._mapper.Mapper._load_table')
    def test_backend_param_passed_to_load(self, mock_load):
        mock_load.return_value = None
        mapper = Mapper()
        Mapper._instance = mapper

        mapper.which_table(
            'uniprot',
            'genesymbol',
            9606,
            backend='biomart',
        )
        mock_load.assert_called_once_with(
            'uniprot',
            'genesymbol',
            9606,
            backend='biomart',
        )

    @patch('omnipath_utils.mapping._mapper.Mapper._load_table')
    def test_backend_forces_reload(self, mock_load):
        """When backend is specified, should reload even if cached."""
        mock_load.return_value = None
        mapper = Mapper()

        # Pre-cache a table
        table = MappingTable(
            data={'P04637': {'TP53'}},
            id_type='uniprot',
            target_id_type='genesymbol',
            ncbi_tax_id=9606,
        )
        mapper.tables[table.key] = table
        Mapper._instance = mapper

        # With backend specified, should force reload
        mapper.which_table(
            'uniprot',
            'genesymbol',
            9606,
            backend='biomart',
        )
        mock_load.assert_called_once()

    @patch('omnipath_utils.mapping._mapper.Mapper._load_table')
    def test_no_backend_uses_cache(self, mock_load):
        """Without backend, cached table should be returned."""
        mapper = Mapper()

        table = MappingTable(
            data={'P04637': {'TP53'}},
            id_type='uniprot',
            target_id_type='genesymbol',
            ncbi_tax_id=9606,
        )
        mapper.tables[table.key] = table
        Mapper._instance = mapper

        result = mapper.which_table('uniprot', 'genesymbol', 9606)
        assert result is table
        mock_load.assert_not_called()


class TestPublicAPIParameters:
    """Test that public API functions accept and pass through new params."""

    def setup_method(self):
        Mapper._instance = None

    @patch('omnipath_utils.mapping._translate.translate_core')
    def test_map_name_passes_raw_and_backend(self, mock_core):
        mock_core.return_value = {'TP53': {'P04637'}}
        from omnipath_utils.mapping import map_name

        map_name('TP53', 'genesymbol', 'uniprot', raw=True, backend='biomart')

        mock_core.assert_called_once()
        call_kwargs = mock_core.call_args
        assert call_kwargs[1]['raw'] is True
        assert call_kwargs[1]['backend'] == 'biomart'

    @patch('omnipath_utils.mapping._translate.translate_core')
    def test_translate_passes_raw_and_backend(self, mock_core):
        mock_core.return_value = {'TP53': {'P04637'}}
        from omnipath_utils.mapping import translate

        translate(
            ['TP53'], 'genesymbol', 'uniprot', raw=True, backend='uniprot'
        )

        mock_core.assert_called_once()
        call_kwargs = mock_core.call_args
        assert call_kwargs[1]['raw'] is True
        assert call_kwargs[1]['backend'] == 'uniprot'

    @patch('omnipath_utils.mapping._translate.translate_core')
    def test_map_names_passes_raw_and_backend(self, mock_core):
        mock_core.return_value = {'TP53': {'P04637'}}
        from omnipath_utils.mapping import map_names

        map_names(
            ['TP53'], 'genesymbol', 'uniprot', raw=True, backend='biomart'
        )

        mock_core.assert_called_once()
        call_kwargs = mock_core.call_args
        assert call_kwargs[1]['raw'] is True
        assert call_kwargs[1]['backend'] == 'biomart'

    @patch('omnipath_utils.mapping._translate.translate_core')
    def test_map_name0_passes_raw_and_backend(self, mock_core):
        mock_core.return_value = {'TP53': {'P04637'}}
        from omnipath_utils.mapping import map_name0

        map_name0('TP53', 'genesymbol', 'uniprot', raw=True, backend='biomart')

        # map_name0 calls map_name which calls translate_core
        mock_core.assert_called()

    @patch('omnipath_utils.mapping._translate.translate_core')
    def test_translate_column_passes_raw_and_backend(self, mock_core):
        import pandas as pd

        mock_core.return_value = {'P04637': {'TP53'}}
        from omnipath_utils.mapping import translate_column

        mapper = Mapper()
        Mapper._instance = mapper

        df = pd.DataFrame({'protein': ['P04637']})
        translate_column(
            df,
            'protein',
            'uniprot',
            'genesymbol',
            raw=True,
            backend='biomart',
        )

        mock_core.assert_called_once()
        call_kwargs = mock_core.call_args
        assert call_kwargs[1]['raw'] is True
        assert call_kwargs[1]['backend'] == 'biomart'


class TestRESTFallbacks:
    """Test the REST API fallback logic."""

    def test_apply_fallbacks_genesymbol_uppercase(self):
        from omnipath_utils.server._routes_mapping import _apply_fallbacks

        mock_session = MagicMock()

        # First call: direct lookup returns nothing for "tp53"
        # The fallback tries uppercase "TP53"
        def mock_translate_ids(session, ids, src, tgt, tax):
            data = {
                'TP53': {'P04637'},
            }
            return ({i: data.get(i, set()) for i in ids if data.get(i)}, {'uniprot'})

        with patch(
            'omnipath_utils.server._routes_mapping.translate_ids',
            side_effect=mock_translate_ids,
        ):
            result = {}
            result = _apply_fallbacks(
                mock_session,
                ['tp53'],
                'genesymbol',
                'uniprot',
                9606,
                result,
            )

        assert result.get('tp53') == {'P04637'}

    def test_apply_fallbacks_chain_via_uniprot(self):
        from omnipath_utils.server._routes_mapping import _apply_fallbacks

        mock_session = MagicMock()
        call_count = {'n': 0}

        def mock_translate_ids(session, ids, src, tgt, tax):
            call_count['n'] += 1
            # genesymbol -> uniprot
            if src == 'entrez' and tgt == 'uniprot':
                return ({'7157': {'P04637'}}, {'uniprot'}) if '7157' in ids else ({}, set())
            # uniprot -> genesymbol
            if src == 'uniprot' and tgt == 'genesymbol':
                return ({'P04637': {'TP53'}}, {'uniprot'}) if 'P04637' in ids else ({}, set())
            return ({}, set())

        with patch(
            'omnipath_utils.server._routes_mapping.translate_ids',
            side_effect=mock_translate_ids,
        ), patch(
            'omnipath_utils.server._routes_mapping.chain_translate',
            return_value=({'7157': {'TP53'}}, {'uniprot+uniprot'}),
        ):
            result = {}
            result = _apply_fallbacks(
                mock_session,
                ['7157'],
                'entrez',
                'genesymbol',
                9606,
                result,
            )

        assert result.get('7157') == {'TP53'}
