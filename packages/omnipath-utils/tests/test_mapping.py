"""Tests for the mapping module."""

import time
from unittest.mock import MagicMock, patch

from omnipath_utils.mapping._table import MappingTable, MappingTableKey
from omnipath_utils.mapping._mapper import Mapper
from omnipath_utils.mapping._cleanup import is_uniprot_ac, uniprot_cleanup
from omnipath_utils.mapping._special import (
    chain_map,
    map_refseq,
    map_genesymbol_fallbacks,
    map_ensembl_strip_version,
)

# ---- MappingTable tests ----


class TestMappingTable:
    def test_create(self, sample_mapping_data):
        table = MappingTable(
            data=sample_mapping_data,
            id_type='uniprot',
            target_id_type='genesymbol',
            ncbi_tax_id=9606,
        )
        assert len(table) == 4

    def test_lookup(self, sample_mapping_data):
        table = MappingTable(
            data=sample_mapping_data,
            id_type='uniprot',
            target_id_type='genesymbol',
            ncbi_tax_id=9606,
        )
        assert table['P04637'] == {'TP53'}
        assert table['NONEXISTENT'] == set()

    def test_contains(self, sample_mapping_data):
        table = MappingTable(
            data=sample_mapping_data,
            id_type='uniprot',
            target_id_type='genesymbol',
            ncbi_tax_id=9606,
        )
        assert 'P04637' in table
        assert 'FAKE' not in table

    def test_key(self, sample_mapping_data):
        table = MappingTable(
            data=sample_mapping_data,
            id_type='uniprot',
            target_id_type='genesymbol',
            ncbi_tax_id=9606,
        )
        assert table.key == MappingTableKey('uniprot', 'genesymbol', 9606)

    def test_expiry(self, sample_mapping_data):
        table = MappingTable(
            data=sample_mapping_data,
            id_type='uniprot',
            target_id_type='genesymbol',
            ncbi_tax_id=9606,
            lifetime=0,
        )
        time.sleep(0.01)
        assert table.expired

    def test_no_expiry(self, sample_mapping_data):
        table = MappingTable(
            data=sample_mapping_data,
            id_type='uniprot',
            target_id_type='genesymbol',
            ncbi_tax_id=9606,
            lifetime=3600,
        )
        assert not table.expired

    def test_repr(self, sample_mapping_data):
        table = MappingTable(
            data=sample_mapping_data,
            id_type='uniprot',
            target_id_type='genesymbol',
            ncbi_tax_id=9606,
        )
        r = repr(table)
        assert 'uniprot' in r
        assert 'genesymbol' in r
        assert '9606' in r

    def test_items_keys_values(self, sample_mapping_data):
        table = MappingTable(
            data=sample_mapping_data,
            id_type='uniprot',
            target_id_type='genesymbol',
            ncbi_tax_id=9606,
        )
        assert set(table.keys()) == {'P04637', 'P00533', 'P38398', 'Q13315'}
        assert {'TP53'} in list(table.values())
        assert len(list(table.items())) == 4

    def test_last_used_updates_on_access(self, sample_mapping_data):
        table = MappingTable(
            data=sample_mapping_data,
            id_type='uniprot',
            target_id_type='genesymbol',
            ncbi_tax_id=9606,
            lifetime=3600,
        )
        before = table._last_used
        time.sleep(0.01)
        _ = table['P04637']
        assert table._last_used > before

    def test_last_used_updates_on_contains(self, sample_mapping_data):
        table = MappingTable(
            data=sample_mapping_data,
            id_type='uniprot',
            target_id_type='genesymbol',
            ncbi_tax_id=9606,
            lifetime=3600,
        )
        before = table._last_used
        time.sleep(0.01)
        _ = 'P04637' in table
        assert table._last_used > before


# ---- MappingTableKey tests ----


class TestMappingTableKey:
    def test_named_tuple(self):
        key = MappingTableKey('uniprot', 'genesymbol', 9606)
        assert key.id_type == 'uniprot'
        assert key.target_id_type == 'genesymbol'
        assert key.ncbi_tax_id == 9606

    def test_equality(self):
        k1 = MappingTableKey('uniprot', 'genesymbol', 9606)
        k2 = MappingTableKey('uniprot', 'genesymbol', 9606)
        assert k1 == k2

    def test_inequality(self):
        k1 = MappingTableKey('uniprot', 'genesymbol', 9606)
        k2 = MappingTableKey('uniprot', 'genesymbol', 10090)
        assert k1 != k2

    def test_hashable(self):
        k = MappingTableKey('uniprot', 'genesymbol', 9606)
        d = {k: 'value'}
        assert d[k] == 'value'


# ---- Cleanup tests ----


class TestUniProtCleanup:
    def test_valid_ac(self):
        assert is_uniprot_ac('P04637')
        assert is_uniprot_ac('Q9Y6K9')
        assert is_uniprot_ac('A0A024R1R8')

    def test_invalid_ac(self):
        assert not is_uniprot_ac('TP53')
        assert not is_uniprot_ac('ENSG00000141510')
        assert not is_uniprot_ac('')
        assert not is_uniprot_ac('12345')

    def test_cleanup_filters_invalid(self):
        mock_mapper = MagicMock()
        mock_mapper._direct_lookup = MagicMock(return_value=set())
        result = uniprot_cleanup(
            {'P04637', 'TP53', 'FAKE'}, 9606, mapper=mock_mapper
        )
        assert 'P04637' in result
        assert 'TP53' not in result
        assert 'FAKE' not in result

    def test_cleanup_empty_input(self):
        result = uniprot_cleanup(set(), 9606)
        assert result == set()


# ---- Special case tests ----


class TestSpecialCases:
    def _make_mock_mapper(self, tables):
        """Create a mock mapper with preloaded tables."""
        mapper = MagicMock()

        def direct_lookup(name, id_type, target_id_type, ncbi_tax_id):
            key = (id_type, target_id_type, ncbi_tax_id)
            data = tables.get(key, {})
            return data.get(name, set())

        mapper._direct_lookup = direct_lookup
        mapper.map_name = MagicMock(return_value=set())
        return mapper

    def test_genesymbol_case_fallback_upper(self):
        tables = {
            ('genesymbol', 'uniprot', 9606): {
                'TP53': {'P04637'},
                'EGFR': {'P00533'},
            },
            ('genesymbol-syn', 'uniprot', 9606): {},
        }
        mapper = self._make_mock_mapper(tables)

        # lowercase should find uppercase via .upper()
        result = map_genesymbol_fallbacks('tp53', 'uniprot', 9606, mapper)
        assert 'P04637' in result

    def test_genesymbol_case_fallback_capitalize(self):
        tables = {
            ('genesymbol', 'uniprot', 10090): {
                'Trp53': {'Q01279'},
            },
            ('genesymbol-syn', 'uniprot', 10090): {},
        }
        mapper = self._make_mock_mapper(tables)

        # all-lower should capitalize and match rodent-style Trp53
        result = map_genesymbol_fallbacks('trp53', 'uniprot', 10090, mapper)
        assert 'Q01279' in result

    def test_genesymbol_synonym_fallback(self):
        tables = {
            ('genesymbol', 'uniprot', 9606): {},
            ('genesymbol-syn', 'uniprot', 9606): {
                'P53': {'P04637'},
            },
        }
        mapper = self._make_mock_mapper(tables)

        result = map_genesymbol_fallbacks('P53', 'uniprot', 9606, mapper)
        assert 'P04637' in result

    def test_genesymbol_strict_mode(self):
        tables = {
            ('genesymbol', 'uniprot', 9606): {},
            ('genesymbol-syn', 'uniprot', 9606): {},
        }
        mapper = self._make_mock_mapper(tables)

        # strict mode should not try appending "1" or truncating
        result = map_genesymbol_fallbacks(
            'SOMEGENE',
            'uniprot',
            9606,
            mapper,
            strict=True,
        )
        assert result == set()

    def test_refseq_version_strip(self):
        tables = {
            ('refseqp', 'uniprot', 9606): {
                'NP_000537': {'P04637'},
            },
        }
        mapper = self._make_mock_mapper(tables)

        # with version should strip and find
        result = map_refseq('NP_000537.3', 'refseqp', 'uniprot', 9606, mapper)
        assert 'P04637' in result

    def test_refseq_exact_match(self):
        tables = {
            ('refseqp', 'uniprot', 9606): {
                'NP_000537.3': {'P04637'},
            },
        }
        mapper = self._make_mock_mapper(tables)

        # exact match with version present
        result = map_refseq(
            'NP_000537.3',
            'refseqp',
            'uniprot',
            9606,
            mapper,
        )
        assert 'P04637' in result

    def test_ensembl_version_strip(self):
        tables = {
            ('ensg', 'genesymbol', 9606): {
                'ENSG00000141510': {'TP53'},
            },
        }
        mapper = self._make_mock_mapper(tables)

        result = map_ensembl_strip_version(
            'ENSG00000141510.12',
            'ensg',
            'genesymbol',
            9606,
            mapper,
        )
        assert 'TP53' in result

    def test_ensembl_no_version_returns_empty(self):
        tables = {
            ('ensg', 'genesymbol', 9606): {},
        }
        mapper = self._make_mock_mapper(tables)

        # No dot in the name -> no stripping, returns empty
        result = map_ensembl_strip_version(
            'ENSG00000141510',
            'ensg',
            'genesymbol',
            9606,
            mapper,
        )
        assert result == set()

    def test_chain_map(self):
        mapper = MagicMock()
        # First call: genesymbol -> uniprot
        # Second call: uniprot -> entrez
        mapper.map_name = MagicMock(
            side_effect=[
                {'P04637'},  # TP53 -> uniprot
                {'7157'},  # P04637 -> entrez
            ],
        )

        result = chain_map('TP53', 'genesymbol', 'entrez', 9606, mapper)
        assert '7157' in result

    def test_chain_map_no_intermediate(self):
        mapper = MagicMock()
        # No intermediate result
        mapper.map_name = MagicMock(return_value=set())

        result = chain_map(
            'NONEXIST',
            'genesymbol',
            'entrez',
            9606,
            mapper,
        )
        assert result == set()


# ---- Mapper tests with mocked backends ----


class TestMapper:
    def setup_method(self):
        """Reset singleton for test isolation."""
        Mapper._instance = None

    def test_mapper_creation(self):
        mapper = Mapper(ncbi_tax_id=9606)
        assert mapper.ncbi_tax_id == 9606
        assert len(mapper.tables) == 0

    def test_default_organism(self):
        mapper = Mapper()
        assert mapper.ncbi_tax_id == 9606

    def test_same_type_returns_self(self):
        mapper = Mapper()
        result = mapper.map_name('P04637', 'uniprot', 'uniprot')
        assert result == {'P04637'}

    def test_empty_name_returns_empty(self):
        mapper = Mapper()
        result = mapper.map_name('', 'uniprot', 'genesymbol')
        assert result == set()

    @patch('omnipath_utils.mapping._mapper.Mapper._load_table')
    def test_map_name_with_preloaded_table(self, mock_load):
        mapper = Mapper()

        # Preload a table
        table = MappingTable(
            data={'P04637': {'TP53'}, 'P00533': {'EGFR'}},
            id_type='uniprot',
            target_id_type='genesymbol',
            ncbi_tax_id=9606,
        )
        mapper.tables[table.key] = table

        result = mapper.map_name('P04637', 'uniprot', 'genesymbol')
        assert result == {'TP53'}

        # Backend should not have been called
        mock_load.assert_not_called()

    @patch('omnipath_utils.mapping._mapper.Mapper._load_table')
    def test_map_names(self, mock_load):
        mapper = Mapper()

        table = MappingTable(
            data={'P04637': {'TP53'}, 'P00533': {'EGFR'}},
            id_type='uniprot',
            target_id_type='genesymbol',
            ncbi_tax_id=9606,
        )
        mapper.tables[table.key] = table

        result = mapper.map_names(
            ['P04637', 'P00533'],
            'uniprot',
            'genesymbol',
        )
        assert 'TP53' in result
        assert 'EGFR' in result

    @patch('omnipath_utils.mapping._mapper.Mapper._load_table')
    def test_translate_batch(self, mock_load):
        mapper = Mapper()

        table = MappingTable(
            data={'P04637': {'TP53'}, 'P00533': {'EGFR'}},
            id_type='uniprot',
            target_id_type='genesymbol',
            ncbi_tax_id=9606,
        )
        mapper.tables[table.key] = table

        result = mapper.translate(
            ['P04637', 'P00533', 'FAKE'],
            'uniprot',
            'genesymbol',
        )
        assert result['P04637'] == {'TP53'}
        assert result['P00533'] == {'EGFR'}
        assert result['FAKE'] == set()

    @patch('omnipath_utils.mapping._mapper.Mapper._load_table')
    def test_map_name0(self, mock_load):
        mapper = Mapper()

        table = MappingTable(
            data={'P04637': {'TP53'}},
            id_type='uniprot',
            target_id_type='genesymbol',
            ncbi_tax_id=9606,
        )
        mapper.tables[table.key] = table

        result = mapper.map_name0('P04637', 'uniprot', 'genesymbol')
        assert result == 'TP53'

    @patch('omnipath_utils.mapping._mapper.Mapper._load_table')
    def test_map_name0_missing(self, mock_load):
        mock_load.return_value = None
        mapper = Mapper()

        result = mapper.map_name0('FAKE', 'uniprot', 'genesymbol')
        assert result is None

    def test_remove_expired(self):
        mapper = Mapper()

        table = MappingTable(
            data={'P04637': {'TP53'}},
            id_type='uniprot',
            target_id_type='genesymbol',
            ncbi_tax_id=9606,
            lifetime=0,
        )
        mapper.tables[table.key] = table

        time.sleep(0.01)
        mapper.remove_expired()
        assert len(mapper.tables) == 0

    def test_remove_expired_keeps_fresh(self):
        mapper = Mapper()

        table = MappingTable(
            data={'P04637': {'TP53'}},
            id_type='uniprot',
            target_id_type='genesymbol',
            ncbi_tax_id=9606,
            lifetime=3600,
        )
        mapper.tables[table.key] = table

        mapper.remove_expired()
        assert len(mapper.tables) == 1

    @patch('omnipath_utils.mapping._mapper.Mapper._load_table')
    def test_which_table_no_load(self, mock_load):
        mapper = Mapper()
        result = mapper.which_table(
            'uniprot',
            'genesymbol',
            9606,
            load=False,
        )
        assert result is None
        mock_load.assert_not_called()

    @patch('omnipath_utils.mapping._mapper.Mapper._load_table')
    def test_translation_table(self, mock_load):
        mapper = Mapper()

        table = MappingTable(
            data={'P04637': {'TP53'}, 'P00533': {'EGFR'}},
            id_type='uniprot',
            target_id_type='genesymbol',
            ncbi_tax_id=9606,
        )
        mapper.tables[table.key] = table

        result = mapper.translation_table('uniprot', 'genesymbol')
        assert isinstance(result, dict)
        assert result['P04637'] == {'TP53'}

    def test_repr(self):
        mapper = Mapper()
        r = repr(mapper)
        assert 'Mapper' in r
        assert '9606' in r

    def test_id_types(self):
        mapper = Mapper()
        types = mapper.id_types()
        assert 'uniprot' in types
        assert 'genesymbol' in types
        assert len(types) > 50


# ---- UniProt Cleanup Pipeline tests ----


class TestUniProtCleanupPipeline:
    """Test the full cleanup pipeline with mock mapper."""

    def _make_mapper_with_tables(self, tables):
        """Create a mock mapper with preloaded lookup tables."""
        mapper = MagicMock()

        def direct_lookup(name, id_type, target_id_type, ncbi_tax_id):
            key = (id_type, target_id_type, ncbi_tax_id)
            data = tables.get(key, {})
            return data.get(name, set())

        mapper._direct_lookup = direct_lookup
        return mapper

    def test_primary_uniprot_translates_secondary(self):
        from omnipath_utils.mapping._cleanup import _primary_uniprot

        mapper = self._make_mapper_with_tables(
            {
                ('uniprot-sec', 'uniprot-pri', 9606): {
                    'Q15086': {'P04637'},  # Q15086 is a secondary AC for TP53
                },
            }
        )
        result = _primary_uniprot({'Q15086', 'P00533'}, 9606, mapper)
        assert 'P04637' in result  # secondary translated
        assert 'P00533' in result  # non-secondary kept

    def test_primary_uniprot_keeps_primary(self):
        from omnipath_utils.mapping._cleanup import _primary_uniprot

        mapper = self._make_mapper_with_tables({})  # no sec->pri table
        result = _primary_uniprot({'P04637'}, 9606, mapper)
        assert result == {'P04637'}

    def test_trembl_to_swissprot(self):
        from omnipath_utils.mapping._cleanup import _trembl_to_swissprot

        mapper = self._make_mapper_with_tables(
            {
                ('trembl', 'genesymbol', 9606): {
                    'A0A024R1R8': {'TP53'},
                },
                ('genesymbol', 'swissprot', 9606): {
                    'TP53': {'P04637'},
                },
            }
        )

        # Mock the reflist to define SwissProt set
        with patch(
            'omnipath_utils.reflists.all_swissprots',
            return_value={'P04637', 'P00533'},
        ):
            result = _trembl_to_swissprot(
                {'A0A024R1R8', 'P00533'}, 9606, mapper
            )

        assert 'P04637' in result  # TrEMBL resolved to SwissProt
        assert 'P00533' in result  # already SwissProt, kept
        assert 'A0A024R1R8' not in result  # TrEMBL replaced

    def test_trembl_kept_when_no_swissprot(self):
        from omnipath_utils.mapping._cleanup import _trembl_to_swissprot

        mapper = self._make_mapper_with_tables(
            {
                ('trembl', 'genesymbol', 9606): {
                    'X12345': {'UNKNOWN_GENE'},
                },
                (
                    'genesymbol',
                    'swissprot',
                    9606,
                ): {},  # no SwissProt for this gene
            }
        )

        with patch(
            'omnipath_utils.reflists.all_swissprots', return_value={'P04637'}
        ):
            result = _trembl_to_swissprot({'X12345'}, 9606, mapper)

        assert 'X12345' in result  # kept because no SwissProt found

    def test_filter_organism(self):
        from omnipath_utils.mapping._cleanup import _filter_organism

        with patch(
            'omnipath_utils.reflists.all_uniprots',
            return_value={'P04637', 'P00533'},
        ):
            result = _filter_organism({'P04637', 'FAKE123'}, 9606)

        assert result == {'P04637'}  # FAKE123 removed

    def test_filter_organism_keeps_all_when_none_match(self):
        from omnipath_utils.mapping._cleanup import _filter_organism

        with patch(
            'omnipath_utils.reflists.all_uniprots', return_value={'OTHER'}
        ):
            result = _filter_organism({'P04637', 'P00533'}, 9606)

        # Nothing matched proteome, but we return originals rather than empty
        assert result == {'P04637', 'P00533'}

    def test_full_cleanup_pipeline(self):
        mapper = self._make_mapper_with_tables(
            {
                ('uniprot-sec', 'uniprot-pri', 9606): {
                    'Q15086': {'P04637'},
                },
            }
        )

        with (
            patch(
                'omnipath_utils.reflists.all_swissprots',
                return_value={'P04637'},
            ),
            patch(
                'omnipath_utils.reflists.all_uniprots',
                return_value={'P04637', 'P00533'},
            ),
        ):
            result = uniprot_cleanup(
                {'Q15086', 'P00533', 'NOT_AC'}, 9606, mapper
            )

        assert 'P04637' in result  # secondary translated
        assert 'P00533' in result  # valid, kept
        assert 'NOT_AC' not in result  # invalid format removed

    def test_cleanup_without_mapper_returns_input(self):
        result = uniprot_cleanup({'P04637', 'FAKE'}, 9606, mapper=None)
        assert result == {'P04637', 'FAKE'}  # returned unchanged


# ---- miRNA fallback tests ----


class TestMirnaFallback:
    def _make_mapper(self, tables):
        mapper = MagicMock()

        def direct_lookup(name, id_type, target_id_type, ncbi_tax_id):
            return tables.get((id_type, target_id_type, ncbi_tax_id), {}).get(
                name, set()
            )

        mapper._direct_lookup = direct_lookup
        return mapper

    def test_mature_to_precursor_fallback(self):
        """mir-mat-name not found, try as mir-name."""
        from omnipath_utils.mapping._special import map_mirna_fallback

        mapper = self._make_mapper(
            {
                # Direct lookup for mir-mat-name fails (empty)
                # But treating it as mir-name succeeds:
                ('mir-name', 'mir-pre', 9606): {
                    'hsa-mir-21': {'MI0000077'},
                },
            }
        )

        result = map_mirna_fallback(
            'hsa-mir-21', 'mir-mat-name', 'mir-pre', 9606, mapper
        )
        assert 'MI0000077' in result

    def test_precursor_to_mature_fallback(self):
        """mir-name not found, try as mir-mat-name."""
        from omnipath_utils.mapping._special import map_mirna_fallback

        mapper = self._make_mapper(
            {
                ('mir-mat-name', 'mirbase', 9606): {
                    'hsa-miR-21-5p': {'MIMAT0000076'},
                },
            }
        )

        result = map_mirna_fallback(
            'hsa-miR-21-5p', 'mir-name', 'mirbase', 9606, mapper
        )
        assert 'MIMAT0000076' in result

    def test_no_mirna_fallback_for_non_mirna(self):
        from omnipath_utils.mapping._special import map_mirna_fallback

        mapper = MagicMock()

        result = map_mirna_fallback(
            'TP53', 'genesymbol', 'uniprot', 9606, mapper
        )
        assert result == set()


# ---- Prefix stripping tests ----


class TestPrefixStripping:
    def test_strip_chebi_prefix(self):
        from omnipath_utils.mapping._special import strip_prefix

        mapper = MagicMock()
        mapper._direct_lookup = MagicMock(return_value={'some_id'})

        result = strip_prefix('CHEBI:12345', 'chebi', 'hmdb', 0, mapper)
        assert result == {'some_id'}
        mapper._direct_lookup.assert_called_with('12345', 'chebi', 'hmdb', 0)

    def test_no_prefix_returns_empty(self):
        from omnipath_utils.mapping._special import strip_prefix

        mapper = MagicMock()

        result = strip_prefix('12345', 'chebi', 'hmdb', 0, mapper)
        assert result == set()


# ---- Backend registration tests ----


class TestBackendRegistry:
    def test_all_backends_register(self):
        from omnipath_utils.mapping.backends import get_backend

        for name in (
            'uniprot',
            'biomart',
            'uploadlists',
            'uniprot_ftp',
            'mirbase',
            'unichem',
            'ramp',
            'hmdb',
        ):
            b = get_backend(name)
            assert b is not None, f'Backend {name} not registered'

    def test_unknown_backend_returns_none(self):
        from omnipath_utils.mapping.backends import get_backend

        assert get_backend('nonexistent') is None

    def test_backend_has_name(self):
        from omnipath_utils.mapping.backends import get_backend

        b = get_backend('uniprot')
        assert b.name == 'uniprot'
        b = get_backend('biomart')
        assert b.name == 'biomart'

    def test_backend_has_yaml_key(self):
        from omnipath_utils.mapping.backends import get_backend

        b = get_backend('uniprot')
        assert b.yaml_key == 'uniprot'
        b = get_backend('biomart')
        assert b.yaml_key == 'ensembl'


# ---- Mapper integration tests (with full fallback chain) ----


class TestMapperIntegration:
    """Test the full map_name pipeline with preloaded tables."""

    def setup_method(self):
        Mapper._instance = None

    @patch('omnipath_utils.mapping._mapper.Mapper._load_table')
    def test_uniprot_target_triggers_cleanup(self, mock_load):
        mapper = Mapper()
        mock_load.return_value = None

        # Preload genesymbol -> uniprot table with a TrEMBL and a SwissProt
        table = MappingTable(
            data={'TP53': {'P04637', 'A0A024R1R8'}},
            id_type='genesymbol',
            target_id_type='uniprot',
            ncbi_tax_id=9606,
        )
        mapper.tables[table.key] = table

        # Without cleanup support (no sec/pri tables, no reflist),
        # cleanup returns the input filtered by format only
        result = mapper.map_name('TP53', 'genesymbol', 'uniprot')
        assert 'P04637' in result
        # Both should be present since there is no reflist to filter TrEMBL
        assert 'A0A024R1R8' in result

    @patch('omnipath_utils.mapping._mapper.Mapper._load_table')
    def test_case_fallback_not_applied_to_uniprot_source(self, mock_load):
        mapper = Mapper()
        mock_load.return_value = None

        # map_name with uniprot source should NOT try uppercase
        result = mapper.map_name('p04637', 'uniprot', 'genesymbol')
        assert result == set()  # no fallback, just empty


# ---- DataFrame translation tests ----


class TestTranslateColumn:
    """Tests for translate_column."""

    def setup_method(self):
        """Reset singleton for test isolation."""
        Mapper._instance = None

    def test_basic_translation(self):
        import pandas as pd

        from omnipath_utils.mapping import translate_column

        mapper = Mapper()
        table = MappingTable(
            data={'P04637': {'TP53'}, 'P00533': {'EGFR'}},
            id_type='uniprot',
            target_id_type='genesymbol',
            ncbi_tax_id=9606,
        )
        mapper.tables[table.key] = table
        Mapper._instance = mapper

        df = pd.DataFrame({'protein': ['P04637', 'P00533', 'FAKE']})
        result = translate_column(df, 'protein', 'uniprot', 'genesymbol')

        assert 'genesymbol' in result.columns
        assert (
            result.loc[result['protein'] == 'P04637', 'genesymbol'].iloc[0]
            == 'TP53'
        )
        assert (
            result.loc[result['protein'] == 'P00533', 'genesymbol'].iloc[0]
            == 'EGFR'
        )

    def test_expand_one_to_many(self):
        import pandas as pd

        from omnipath_utils.mapping import translate_column

        mapper = Mapper()
        table = MappingTable(
            data={'TP53': {'P04637', 'A0A024R1R8'}},
            id_type='genesymbol',
            target_id_type='uniprot',
            ncbi_tax_id=9606,
        )
        mapper.tables[table.key] = table
        Mapper._instance = mapper

        df = pd.DataFrame({'gene': ['TP53']})
        result = translate_column(
            df,
            'gene',
            'genesymbol',
            'uniprot',
            expand=True,
        )

        # Should have 2 rows (one per UniProt)
        assert len(result) == 2
        assert set(result['uniprot']) == {'P04637', 'A0A024R1R8'}

    def test_no_expand(self):
        import pandas as pd

        from omnipath_utils.mapping import translate_column

        mapper = Mapper()
        table = MappingTable(
            data={'TP53': {'P04637', 'A0A024R1R8'}},
            id_type='genesymbol',
            target_id_type='uniprot',
            ncbi_tax_id=9606,
        )
        mapper.tables[table.key] = table
        Mapper._instance = mapper

        df = pd.DataFrame({'gene': ['TP53']})
        result = translate_column(
            df,
            'gene',
            'genesymbol',
            'uniprot',
            expand=False,
        )

        # Should have 1 row (picked one)
        assert len(result) == 1
        assert result['uniprot'].iloc[0] in {'P04637', 'A0A024R1R8'}

    def test_drop_untranslated(self):
        import pandas as pd

        from omnipath_utils.mapping import translate_column

        mapper = Mapper()
        table = MappingTable(
            data={'P04637': {'TP53'}},
            id_type='uniprot',
            target_id_type='genesymbol',
            ncbi_tax_id=9606,
        )
        mapper.tables[table.key] = table
        Mapper._instance = mapper

        df = pd.DataFrame({'protein': ['P04637', 'FAKE']})
        result = translate_column(
            df,
            'protein',
            'uniprot',
            'genesymbol',
            keep_untranslated=False,
        )

        assert len(result) == 1
        assert result['protein'].iloc[0] == 'P04637'

    def test_keep_untranslated(self):
        import pandas as pd

        from omnipath_utils.mapping import translate_column

        mapper = Mapper()
        table = MappingTable(
            data={'P04637': {'TP53'}},
            id_type='uniprot',
            target_id_type='genesymbol',
            ncbi_tax_id=9606,
        )
        mapper.tables[table.key] = table
        Mapper._instance = mapper

        df = pd.DataFrame({'protein': ['P04637', 'FAKE']})
        result = translate_column(
            df,
            'protein',
            'uniprot',
            'genesymbol',
            keep_untranslated=True,
        )

        assert len(result) == 2
        assert pd.isna(
            result.loc[result['protein'] == 'FAKE', 'genesymbol'].iloc[0]
        )

    def test_custom_new_column(self):
        import pandas as pd

        from omnipath_utils.mapping import translate_column

        mapper = Mapper()
        table = MappingTable(
            data={'P04637': {'TP53'}},
            id_type='uniprot',
            target_id_type='genesymbol',
            ncbi_tax_id=9606,
        )
        mapper.tables[table.key] = table
        Mapper._instance = mapper

        df = pd.DataFrame({'protein': ['P04637']})
        result = translate_column(
            df,
            'protein',
            'uniprot',
            'genesymbol',
            new_column='gene_name',
        )

        assert 'gene_name' in result.columns
        assert 'genesymbol' not in result.columns

    def test_does_not_modify_input(self):
        import pandas as pd

        from omnipath_utils.mapping import translate_column

        mapper = Mapper()
        table = MappingTable(
            data={'P04637': {'TP53'}},
            id_type='uniprot',
            target_id_type='genesymbol',
            ncbi_tax_id=9606,
        )
        mapper.tables[table.key] = table
        Mapper._instance = mapper

        df = pd.DataFrame({'protein': ['P04637']})
        original_cols = list(df.columns)
        _ = translate_column(df, 'protein', 'uniprot', 'genesymbol')

        assert list(df.columns) == original_cols


class TestTranslateColumns:
    """Tests for translate_columns."""

    def setup_method(self):
        """Reset singleton for test isolation."""
        Mapper._instance = None

    def test_multi_translation(self):
        import pandas as pd

        from omnipath_utils.mapping import translate_columns

        mapper = Mapper()

        t1 = MappingTable(
            data={'P04637': {'TP53'}},
            id_type='uniprot',
            target_id_type='genesymbol',
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

        df = pd.DataFrame({'protein': ['P04637']})
        result = translate_columns(
            df,
            ('protein', 'uniprot', 'genesymbol'),
            ('protein', 'uniprot', 'entrez', 'entrez_id'),
            expand=False,
        )

        assert 'genesymbol' in result.columns
        assert 'entrez_id' in result.columns
        assert result['genesymbol'].iloc[0] == 'TP53'
        assert result['entrez_id'].iloc[0] == '7157'

    def test_chained_expand(self):
        import pandas as pd

        from omnipath_utils.mapping import translate_columns

        mapper = Mapper()

        t1 = MappingTable(
            data={'P04637': {'TP53'}},
            id_type='uniprot',
            target_id_type='genesymbol',
            ncbi_tax_id=9606,
        )
        t2 = MappingTable(
            data={'P04637': {'7157', '99999'}},
            id_type='uniprot',
            target_id_type='entrez',
            ncbi_tax_id=9606,
        )
        mapper.tables[t1.key] = t1
        mapper.tables[t2.key] = t2
        Mapper._instance = mapper

        df = pd.DataFrame({'protein': ['P04637']})
        result = translate_columns(
            df,
            ('protein', 'uniprot', 'genesymbol'),
            ('protein', 'uniprot', 'entrez'),
            expand=True,
        )

        # 1 genesymbol * 2 entrez = 2 rows
        assert len(result) == 2
        assert set(result['entrez']) == {'7157', '99999'}


class TestTranslateColumnPolars:
    """Tests for translate_column with polars backend."""

    def setup_method(self):
        """Reset singleton for test isolation."""
        Mapper._instance = None

    def test_basic_with_polars(self):
        pl = __import__('pytest').importorskip('polars')
        from omnipath_utils.mapping import translate_column

        mapper = Mapper()
        table = MappingTable(
            data={'P04637': {'TP53'}, 'P00533': {'EGFR'}},
            id_type='uniprot',
            target_id_type='genesymbol',
            ncbi_tax_id=9606,
        )
        mapper.tables[table.key] = table
        Mapper._instance = mapper

        df = pl.DataFrame({'protein': ['P04637', 'P00533', 'FAKE']})
        result = translate_column(df, 'protein', 'uniprot', 'genesymbol')

        assert isinstance(result, pl.DataFrame)
        assert 'genesymbol' in result.columns

    def test_expand_with_polars(self):
        pl = __import__('pytest').importorskip('polars')
        from omnipath_utils.mapping import translate_column

        mapper = Mapper()
        table = MappingTable(
            data={'TP53': {'P04637', 'A0A024R1R8'}},
            id_type='genesymbol',
            target_id_type='uniprot',
            ncbi_tax_id=9606,
        )
        mapper.tables[table.key] = table
        Mapper._instance = mapper

        df = pl.DataFrame({'gene': ['TP53']})
        result = translate_column(
            df, 'gene', 'genesymbol', 'uniprot', expand=True
        )

        assert isinstance(result, pl.DataFrame)
        assert len(result) == 2
        assert set(result['uniprot'].to_list()) == {'P04637', 'A0A024R1R8'}
