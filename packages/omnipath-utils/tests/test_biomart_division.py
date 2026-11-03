"""Division-aware BioMart endpoint routing (007 US1 / T017)."""

from __future__ import annotations

import pytest

from omnipath_utils.mapping.backends._biomart import (
    BIOMART_URL,
    ENSEMBL_GENOMES_DIVISIONS,
    BioMartBackend,
)


def test_vertebrate_uses_main_host():
    t = BioMartBackend._target(9606)  # human
    assert t is not None
    assert t.host == BIOMART_URL
    assert t.schema == 'default'
    assert t.dataset.endswith('_gene_ensembl')
    assert t.division is None


def test_classic_model_organisms_use_main_host():
    # fly / worm / yeast live on the MAIN Ensembl host, not a genomes division.
    for taxon in (7227, 6239, 4932):
        t = BioMartBackend._target(taxon)
        assert t is not None, taxon
        assert t.host == BIOMART_URL, taxon
        assert t.dataset.endswith('_gene_ensembl'), taxon
        assert t.division is None, taxon


def test_plant_routes_to_genomes_division():
    t = BioMartBackend._target(3702)  # Arabidopsis thaliana
    assert t is not None
    assert t.division == 'plants'
    host, schema = ENSEMBL_GENOMES_DIVISIONS['plants']
    assert t.host == host
    assert t.schema == schema
    assert t.dataset.endswith('_eg_gene')


def test_unknown_organism_returns_none():
    assert BioMartBackend._target(0) is None
