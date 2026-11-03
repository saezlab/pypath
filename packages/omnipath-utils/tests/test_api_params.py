"""Tests for unified API features: raw mode, backend selection, swissprot/trembl targets."""

from unittest.mock import patch

import pytest

from omnipath_utils.mapping._table import MappingTable
from omnipath_utils.mapping._mapper import Mapper

# ---- Raw mode tests ----


class TestRawMode:
    """Test that raw=True bypasses all special-case handling."""

    def setup_method(self):
        Mapper._instance = None

    @patch('omnipath_utils.mapping._mapper.Mapper._load_table')
    def test_raw_skips_gene_symbol_case_fallback(self, mock_load):
        """In raw mode, 'tp53' (lowercase) should NOT find TP53."""
        from omnipath_utils.mapping import map_name

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

        # Normal mode: lowercase should find via UPPER fallback
        assert 'P04637' in map_name('tp53', 'genesymbol', 'uniprot')

        # Raw mode: no fallback, so lowercase misses
        assert map_name('tp53', 'genesymbol', 'uniprot', raw=True) == set()

    @patch('omnipath_utils.mapping._mapper.Mapper._load_table')
    def test_raw_skips_chain_translation(self, mock_load):
        """In raw mode, chain translation should not happen."""
        from omnipath_utils.mapping import map_name

        mock_load.return_value = None
        mapper = Mapper()
        # genesymbol->uniprot and uniprot->entrez exist, but genesymbol->entrez doesn't
        t1 = MappingTable(
            data={'TP53': {'P04637'}},
            id_type='genesymbol',
            target_id_type='uniprot',
            ncbi_tax_id=9606,
        )
        t2 = MappingTable(
            data={'P04637': {'7157'}},
            id_type='uniprot',
            target_id_type='entrez',
            ncbi_tax_id=9606,
        )
        mapper.tables[t1.key] = t1
        mapper.tables[t2.key] = t2
        Mapper._instance = mapper

        # Normal mode: chains via uniprot
        assert '7157' in map_name('TP53', 'genesymbol', 'entrez')

        # Raw mode: no chain, no direct table -> empty
        assert map_name('TP53', 'genesymbol', 'entrez', raw=True) == set()

    @patch('omnipath_utils.mapping._mapper.Mapper._load_table')
    def test_raw_skips_uniprot_cleanup(self, mock_load):
        """In raw mode, UniProt cleanup should not run."""
        from omnipath_utils.mapping import map_name

        mock_load.return_value = None
        mapper = Mapper()
        # Table returns a non-AC string -- cleanup would filter it out
        table = MappingTable(
            data={'GENE1': {'NOT_AN_AC', 'P04637'}},
            id_type='genesymbol',
            target_id_type='uniprot',
            ncbi_tax_id=9606,
        )
        mapper.tables[table.key] = table
        Mapper._instance = mapper

        # Raw mode: both returned (no format filtering)
        result = map_name('GENE1', 'genesymbol', 'uniprot', raw=True)
        assert 'NOT_AN_AC' in result
        assert 'P04637' in result

    @patch('omnipath_utils.mapping._mapper.Mapper._load_table')
    def test_raw_direct_hit(self, mock_load):
        """In raw mode, direct hits still work."""
        from omnipath_utils.mapping import map_name

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

        assert map_name('TP53', 'genesymbol', 'uniprot', raw=True) == {'P04637'}

    @patch('omnipath_utils.mapping._mapper.Mapper._load_table')
    def test_raw_returns_empty_for_missing_table(self, mock_load):
        """In raw mode with no matching table, returns empty set."""
        from omnipath_utils.mapping import map_name

        mock_load.return_value = None
        mapper = Mapper()
        Mapper._instance = mapper

        result = map_name('TP53', 'genesymbol', 'uniprot', raw=True)
        assert result == set()


# ---- Backend selection tests ----


class TestBackendSelection:
    """Test explicit backend selection."""

    def setup_method(self):
        Mapper._instance = None

    @patch('omnipath_utils.mapping._mapper.Mapper._load_table')
    def test_backend_passed_to_load(self, mock_load):
        from omnipath_utils.mapping import map_name

        mock_load.return_value = None
        mapper = Mapper()
        Mapper._instance = mapper

        map_name('TP53', 'genesymbol', 'uniprot', backend='biomart')

        # _load_table should have been called with backend='biomart'
        calls = mock_load.call_args_list
        assert any(
            c.kwargs.get('backend') == 'biomart'
            or (len(c.args) > 3 and c.args[3] == 'biomart')
            for c in calls
        )

    @patch('omnipath_utils.mapping._mapper.Mapper._load_table')
    def test_backend_in_translate(self, mock_load):
        from omnipath_utils.mapping import translate

        mock_load.return_value = None
        mapper = Mapper()
        Mapper._instance = mapper

        # No tables preloaded -- with raw=True and no table, should return empty
        result = translate(
            ['TP53'],
            'genesymbol',
            'uniprot',
            raw=True,
            backend='nonexistent',
        )
        assert result == {'TP53': set()}

    @patch('omnipath_utils.mapping._mapper.Mapper._load_table')
    def test_backend_forces_reload(self, mock_load):
        """When backend is specified, should bypass cached table and reload."""
        mock_load.return_value = None
        mapper = Mapper()

        # Pre-cache a table (auto-loaded, no specific backend)
        table = MappingTable(
            data={'P04637': {'TP53'}},
            id_type='uniprot',
            target_id_type='genesymbol',
            ncbi_tax_id=9606,
        )
        mapper.tables[table.key] = table
        Mapper._instance = mapper

        # With backend specified, should call _load_table even though cached
        mapper.which_table(
            'uniprot',
            'genesymbol',
            9606,
            backend='biomart',
        )
        mock_load.assert_called_once()


# ---- SwissProt / TrEMBL target types ----


class TestSwissProtTremblTargets:
    """Test that swissprot and trembl work as target types."""

    def setup_method(self):
        Mapper._instance = None

    @patch('omnipath_utils.mapping._mapper.Mapper._load_table')
    def test_genesymbol_to_swissprot(self, mock_load):
        from omnipath_utils.mapping import map_name

        mock_load.return_value = None
        mapper = Mapper()
        table = MappingTable(
            data={'TP53': {'P04637'}},
            id_type='genesymbol',
            target_id_type='swissprot',
            ncbi_tax_id=9606,
        )
        mapper.tables[table.key] = table
        Mapper._instance = mapper

        result = map_name('TP53', 'genesymbol', 'swissprot')
        assert result == {'P04637'}

    @patch('omnipath_utils.mapping._mapper.Mapper._load_table')
    def test_genesymbol_to_trembl(self, mock_load):
        from omnipath_utils.mapping import map_name

        mock_load.return_value = None
        mapper = Mapper()
        table = MappingTable(
            data={'TP53': {'A0A024R1R8'}},
            id_type='genesymbol',
            target_id_type='trembl',
            ncbi_tax_id=9606,
        )
        mapper.tables[table.key] = table
        Mapper._instance = mapper

        result = map_name('TP53', 'genesymbol', 'trembl')
        assert result == {'A0A024R1R8'}


# ---- translate_column with new params ----


class TestTranslateColumnParams:
    def setup_method(self):
        Mapper._instance = None

    @patch('omnipath_utils.mapping._mapper.Mapper._load_table')
    def test_translate_column_raw(self, mock_load):
        pd = pytest.importorskip('pandas', reason='pandas not available')
        from omnipath_utils.mapping import translate_column

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

        df = pd.DataFrame({'gene': ['TP53', 'tp53']})

        # Raw mode: tp53 (lowercase) should not match
        result = translate_column(
            df,
            'gene',
            'genesymbol',
            'uniprot',
            raw=True,
            expand=False,
        )
        assert (
            result.loc[result['gene'] == 'TP53', 'uniprot'].iloc[0] == 'P04637'
        )
        assert pd.isna(result.loc[result['gene'] == 'tp53', 'uniprot'].iloc[0])

    @patch('omnipath_utils.mapping._mapper.Mapper._load_table')
    def test_translate_column_backend(self, mock_load):
        pd = pytest.importorskip('pandas', reason='pandas not available')
        from omnipath_utils.mapping import translate_column

        mock_load.return_value = None
        mapper = Mapper()
        Mapper._instance = mapper

        # With nonexistent backend and raw=True, should get None for all
        df = pd.DataFrame({'gene': ['TP53']})
        result = translate_column(
            df,
            'gene',
            'genesymbol',
            'uniprot',
            raw=True,
            backend='nonexistent',
            expand=False,
        )
        assert pd.isna(result['uniprot'].iloc[0])
