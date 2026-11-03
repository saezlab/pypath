"""Tests for the ID type registry."""

import pytest

from omnipath_utils.mapping._id_types import IdTypeRegistry


@pytest.fixture
def registry():
    return IdTypeRegistry.get()


class TestIdTypeRegistry:
    def test_loads_types(self, registry):
        assert len(registry) > 50

    def test_resolve_canonical(self, registry):
        assert registry.resolve('uniprot') == 'uniprot'
        assert registry.resolve('genesymbol') == 'genesymbol'

    def test_resolve_alias(self, registry):
        assert registry.resolve('gene_symbol') == 'genesymbol'
        assert registry.resolve('uniprotkb') == 'uniprot'

    def test_resolve_dash_underscore(self, registry):
        # genesymbol-syn and genesymbol_syn should both resolve
        result = registry.resolve('genesymbol-syn')
        assert result is not None

    def test_resolve_underscore_variant(self, registry):
        result = registry.resolve('genesymbol_syn')
        assert result is not None

    def test_resolve_unknown(self, registry):
        assert registry.resolve('completely_fake_type') is None

    def test_entity_type(self, registry):
        assert registry.entity_type('uniprot') == 'protein'
        assert registry.entity_type('hmdb') == 'small_molecule'
        assert registry.entity_type('ensg') == 'gene'

    def test_curie_prefix(self, registry):
        assert registry.curie_prefix('uniprot') == 'uniprot'
        assert registry.curie_prefix('chebi') == 'chebi'
        assert registry.curie_prefix('genesymbol') == 'hgnc.symbol'

    def test_backend_column(self, registry):
        assert registry.backend_column('uniprot', 'uniprot') == 'accession'
        assert (
            registry.backend_column('genesymbol', 'uniprot') == 'gene_primary'
        )
        assert (
            registry.backend_column('genesymbol', 'ensembl')
            == 'external_gene_name'
        )

    def test_by_entity_type(self, registry):
        proteins = registry.by_entity_type('protein')
        assert 'uniprot' in proteins
        assert 'swissprot' in proteins
        assert len(proteins) >= 5

    def test_by_backend(self, registry):
        uniprot_types = registry.by_backend('uniprot')
        assert 'genesymbol' in uniprot_types
        assert uniprot_types['genesymbol'] == 'gene_primary'

    def test_all_names(self, registry):
        names = registry.all_names()
        assert 'uniprot' in names
        assert 'genesymbol' in names
        assert len(names) > 50

    def test_contains(self, registry):
        assert 'uniprot' in registry
        assert 'fake_type' not in registry


    def test_chembl_type(self, registry):
        assert registry.resolve('chembl') == 'chembl'
        assert registry.resolve('chembl_id') == 'chembl'
        assert registry.resolve('chembl_compound') == 'chembl'
        assert registry.entity_type('chembl') == 'small_molecule'
        assert registry.curie_prefix('chembl') == 'chembl.compound'
        assert registry.backend_column('chembl', 'unichem') == 'chembl'

    def test_singleton(self):
        r1 = IdTypeRegistry.get()
        r2 = IdTypeRegistry.get()
        assert r1 is r2


class TestAutoDiscoveryTypes:
    """New UniChem and RaMP types are registered in id_types.yaml."""

    def test_unichem_types_registered(self, registry):
        """UniChem sources should all be in the registry."""
        unichem_types = [
            'chembl',
            'drugbank',
            'chebi',
            'hmdb',
            'pubchem',
            'guide_to_pharmacology',
            'swisslipids',
            'bindingdb',
            'drugcentral',
            'surechembl',
            'molport',
            'nmrshiftdb',
            'fda_srs',
            'probes_drugs',
            'csd',
            'rscb_pdb',
            'pdbe',
        ]

        for t in unichem_types:
            assert t in registry, f'{t} not in registry'

    def test_ramp_types_in_registry(self, registry):
        """RaMP compound types should be in the registry."""
        ramp_types = [
            'wikidata',
            'refmet',
            'chemspider',
            'lipidbank',
            'plantfa',
            'kegg_glycan',
            'polymer',
        ]

        for t in ramp_types:
            assert t in registry, f'{t} not in registry'

    def test_unichem_entity_types(self, registry):
        """All UniChem types should be small_molecule."""
        for t in [
            'guide_to_pharmacology',
            'surechembl',
            'bindingdb',
            'drugcentral',
        ]:
            assert registry.entity_type(t) == 'small_molecule'

    def test_unichem_aliases(self, registry):
        """UniChem type aliases should resolve correctly."""
        assert registry.resolve('gtopdb') == 'guide_to_pharmacology'
        assert registry.resolve('gtop') == 'guide_to_pharmacology'
        assert registry.resolve('surechembl_id') == 'surechembl'
        assert registry.resolve('bindingdb_id') == 'bindingdb'
        assert registry.resolve('drugcentral_id') == 'drugcentral'
        assert registry.resolve('probes_and_drugs') == 'probes_drugs'

    def test_refmet_type(self, registry):
        """RefMet should be registered from RaMP."""
        assert 'refmet' in registry
        assert registry.entity_type('refmet') == 'small_molecule'
        assert registry.resolve('refmet_id') == 'refmet'


class TestAutoDiscoveryBuild:
    """Test the _build.py auto-discovery helper methods."""

    def test_unichem_name_map(self):
        """Verify the UniChem name mapping table."""
        from unittest.mock import patch

        with (
            patch('omnipath_utils.db._build.get_engine'),
            patch('omnipath_utils.db._build.ensure_schema'),
        ):
            from omnipath_utils.db._build import DatabaseBuilder

            builder = DatabaseBuilder(db_url='postgresql://test/db')

        assert builder._UNICHEM_NAME_MAP['lipid_maps'] == 'lipidmaps'
        assert builder._UNICHEM_NAME_MAP['probes&drugs'] == 'probes_drugs'

    def test_unichem_canonical(self):
        """Verify canonical name normalisation for UniChem labels."""
        from unittest.mock import patch

        with (
            patch('omnipath_utils.db._build.get_engine'),
            patch('omnipath_utils.db._build.ensure_schema'),
        ):
            from omnipath_utils.db._build import DatabaseBuilder

            builder = DatabaseBuilder(db_url='postgresql://test/db')

        assert builder._unichem_canonical('ChEMBL') == 'chembl'
        assert builder._unichem_canonical('DrugBank') == 'drugbank'
        assert (
            builder._unichem_canonical('Guide to Pharmacology')
            == 'guide_to_pharmacology'
        )
        assert builder._unichem_canonical('LIPID MAPS\u00ae') == 'lipidmaps'
        assert builder._unichem_canonical('Probes&Drugs') == 'probes_drugs'
        assert builder._unichem_canonical('') is None

    def test_ramp_name_map(self):
        """Verify the RaMP name mapping table."""
        from unittest.mock import patch

        with (
            patch('omnipath_utils.db._build.get_engine'),
            patch('omnipath_utils.db._build.ensure_schema'),
        ):
            from omnipath_utils.db._build import DatabaseBuilder

            builder = DatabaseBuilder(db_url='postgresql://test/db')

        assert builder._RAMP_NAME_MAP['CAS'] == 'cas'
        assert builder._RAMP_NAME_MAP['LIPIDMAPS'] == 'lipidmaps'
        assert builder._RAMP_NAME_MAP['rhea-comp'] == 'rhea'

    def test_ramp_canonical(self):
        """Verify canonical name normalisation for RaMP types."""
        from unittest.mock import patch

        with (
            patch('omnipath_utils.db._build.get_engine'),
            patch('omnipath_utils.db._build.ensure_schema'),
        ):
            from omnipath_utils.db._build import DatabaseBuilder

            builder = DatabaseBuilder(db_url='postgresql://test/db')

        assert builder._ramp_canonical('chebi') == 'chebi'
        assert builder._ramp_canonical('CAS') == 'cas'
        assert builder._ramp_canonical('LIPIDMAPS') == 'lipidmaps'
        assert builder._ramp_canonical('rhea-comp') == 'rhea'
        assert builder._ramp_canonical('hmdb') == 'hmdb'
