"""Tests for the taxonomy module."""

from omnipath_utils.taxonomy import (
    all_organisms,
    ensure_oma_code,
    ensure_kegg_code,
    ensure_latin_name,
    ensure_common_name,
    ensure_ncbi_tax_id,
    ensure_ensembl_name,
    ensure_mirbase_name,
)
from omnipath_utils.taxonomy._taxonomy import TaxonomyManager


class TestEnsureNcbiTaxId:
    """Test ensure_ncbi_tax_id with various input formats."""

    def test_integer_input(self):
        assert ensure_ncbi_tax_id(9606) == 9606

    def test_string_integer(self):
        assert ensure_ncbi_tax_id('9606') == 9606

    def test_common_name(self):
        assert ensure_ncbi_tax_id('human') == 9606

    def test_common_name_case_insensitive(self):
        assert ensure_ncbi_tax_id('Human') == 9606
        assert ensure_ncbi_tax_id('HUMAN') == 9606

    def test_latin_name(self):
        assert ensure_ncbi_tax_id('Homo sapiens') == 9606

    def test_ensembl_name(self):
        assert ensure_ncbi_tax_id('hsapiens') == 9606

    def test_kegg_code(self):
        assert ensure_ncbi_tax_id('hsa') == 9606

    def test_mirbase_code(self):
        # hsa is both kegg and mirbase for human
        assert ensure_ncbi_tax_id('hsa') == 9606

    def test_oma_code(self):
        # HUMAN lowered to human matches common_name; the lookup is
        # case-insensitive and oma_code "HUMAN" is indexed lower-cased
        assert ensure_ncbi_tax_id('HUMAN') == 9606

    def test_mouse(self):
        assert ensure_ncbi_tax_id('mouse') == 10090
        assert ensure_ncbi_tax_id(10090) == 10090
        assert ensure_ncbi_tax_id('mmusculus') == 10090

    def test_rat(self):
        assert ensure_ncbi_tax_id('rat') == 10116

    def test_fly(self):
        assert ensure_ncbi_tax_id('fruit fly') == 7227
        assert ensure_ncbi_tax_id('dme') == 7227

    def test_unknown_returns_none(self):
        assert ensure_ncbi_tax_id('alien') is None
        assert ensure_ncbi_tax_id(999999999) is None


class TestEnsureNames:
    """Test name conversion functions."""

    def test_common_name(self):
        assert ensure_common_name(9606) == 'human'
        assert ensure_common_name(10090) == 'mouse'

    def test_latin_name(self):
        assert ensure_latin_name(9606) == 'Homo sapiens'
        assert ensure_latin_name(10090) == 'Mus musculus'

    def test_ensembl_name(self):
        assert ensure_ensembl_name(9606) == 'hsapiens'

    def test_kegg_code(self):
        assert ensure_kegg_code(9606) == 'hsa'
        assert ensure_kegg_code(10090) == 'mmu'
        assert ensure_kegg_code(10116) == 'rno'

    def test_mirbase_name(self):
        assert ensure_mirbase_name(9606) == 'hsa'
        assert ensure_mirbase_name(10090) == 'mmu'

    def test_oma_code(self):
        assert ensure_oma_code(9606) == 'HUMAN'
        assert ensure_oma_code(7227) == 'DROME'

    def test_unknown_returns_none(self):
        assert ensure_common_name(999999999) is None


class TestAllOrganisms:
    def test_returns_dict(self):
        orgs = all_organisms()
        assert isinstance(orgs, dict)
        assert 9606 in orgs
        assert len(orgs) >= 20

    def test_organism_has_fields(self):
        orgs = all_organisms()
        human = orgs[9606]
        assert human['common_name'] == 'human'
        assert human['latin_name'] == 'Homo sapiens'
        assert human['kegg_code'] == 'hsa'


class TestSingleton:
    def test_singleton(self):
        tm1 = TaxonomyManager.get()
        tm2 = TaxonomyManager.get()
        assert tm1 is tm2

    def test_contains(self):
        tm = TaxonomyManager.get()
        assert 9606 in tm
        assert 'human' in tm
        assert 'alien' not in tm

    def test_len(self):
        tm = TaxonomyManager.get()
        assert len(tm) >= 20


class TestSynonyms:
    """Test synonym-based organism lookups."""

    def test_worm(self):
        assert ensure_ncbi_tax_id('worm') == 6239

    def test_fly(self):
        assert ensure_ncbi_tax_id('fly') == 7227

    def test_fish(self):
        assert ensure_ncbi_tax_id('fish') == 7955

    def test_drosophila(self):
        assert ensure_ncbi_tax_id('drosophila') == 7227

    def test_bovine(self):
        assert ensure_ncbi_tax_id('bovine') == 9913

    def test_fission_yeast(self):
        assert ensure_ncbi_tax_id('fission yeast') == 284812

    def test_ecoli(self):
        assert ensure_ncbi_tax_id('e.coli') == 83333

    def test_parenthetical(self):
        assert ensure_ncbi_tax_id('Caenorhabditis elegans (PRJNA13758)') == 6239

    def test_case_insensitive_synonym(self):
        assert ensure_ncbi_tax_id('DROSOPHILA') == 7227
