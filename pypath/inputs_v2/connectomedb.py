"""
Parse ConnectomeDB2025 data and emit Entity records.

This module converts ConnectomeDB2025 interactions and complexes into Entity 
records using the declarative schema pattern.
"""

from __future__ import annotations

import re

from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    InteractionMetadataCv,
    LicenseCV,
    UpdateCategoryCV,
    ResourceCv,
    InterCellAnnotations,
)
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    Member,
    MembershipBuilder,
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.base import iter_csv


# =============================================================================
# Resource Configuration
# =============================================================================

config = ResourceConfig(
    id=ResourceCv.CONNECTOMEDB,
    name='ConnectomeDB2025',
    url='https://connectomedb.org/',
    license=LicenseCV.CC_BY_NC_4_0,
    update_category=UpdateCategoryCV.IRREGULAR,
    pubmed='41171146',
    primary_category='interactions',
    description=(
        'ConnectomeDB is a comprehensive and ongoing project that provides a'
        'high-quality manually curated database of interacting ligand-receptor'
        'pairs for use in cell-cell communication analysis. First released in'
        '2015 (Ramilowski, et al.), and subsequently updated in 2020 (Hou, et al.)'
        'and 2025 (Liu, Maezono et al.), it aims to enhance the understanding of'
        'cell-cell communication in humans and other mammals, supporting biological'
        'and medical research.'
    ),
)


# =============================================================================
# Download Configurations
# =============================================================================

BASE_URL = 'https://connectomedb.org/downloads/Current-Release/CSV/'

download_interactions = Download(
    url=BASE_URL + 'all_species.csv',
    filename='connectomedb_all_species_interactions.csv',
    subfolder='connectomedb2025',
)


# =============================================================================
# Processing Helpers
# =============================================================================

_symbols_pat = re.compile(r"^([^,(]+)(?:\s*\((.+)\))?")


def _extract_primary_gene(token: str):
    match = _symbols_pat.search(token or '')
    return match.group(1).strip() if match else None


def _extract_gene_alias(token: str):
    match = _symbols_pat.search(token or '')
    return match.group(2).replace(", ", ";") if match and match.group(2) else None


_species_taxon = {
    'human': '9606',
    'mouse': '10090',
    'chimp': '9598',
    'macaque': '9544',
    'marmoset': '9483',
    'rat': '10116',
    'pig': '9823',
    'cow': '9913',
    'dog': '9615',
    'horse': '9796',
    'sheep': '9940',
    'chicken': '9031',
    'frog': '8364',
    'zebrafish': '7955',
}


def _species_to_taxon(species: str) -> str | None:
    return _species_taxon.get((species or '').strip().lower())

_cdb_pat = re.compile(r"^CDB\d{2}:(\d+)", re.IGNORECASE)
_hgnc_pat = re.compile(r"^HGNC:(\d+)$", re.IGNORECASE)

def _extract_cdb(val: str) -> str | None:
    return val if _cdb_pat.match(val) else None


def _extract_hgnc_id(val: str) -> str | None:
    match = _hgnc_pat.match((val or '').strip())
    return match.group(1) if match else None


def _location_terms(location: str | None) -> list[InterCellAnnotations]:
    """Map ConnectomeDB location labels to UniProt keyword CV terms."""
    location_lower = (location or '').lower()
    terms: list[InterCellAnnotations] = []

    if 'secreted' in location_lower:
        terms.append(InterCellAnnotations.SECRETED)
    if 'membrane' in location_lower:
        terms.append(InterCellAnnotations.MEMBRANE)
    if 'cytoplasm' in location_lower:
        terms.append(InterCellAnnotations.CYTOPLASM)

    return terms


# =============================================================================
# Field and Schema Definitions
# =============================================================================

f = FieldConfig(
    extract={
        'primary_gene': _extract_primary_gene,
        'gene_alias': _extract_gene_alias,
        'cdb': _extract_cdb,
        'hgnc_id': _extract_hgnc_id,
    },
    map={
        'species_taxon': _species_to_taxon,
    },
)

# -----------------------------------------------------------------------------
# Interactions Schema
# -----------------------------------------------------------------------------


interactions_schema = EntityBuilder(
    entity_type=EntityTypeCv.INTERACTION,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.CDB, value=f('Interaction ID', extract='cdb')),
        CV(term=IdentifierNamespaceCv.NAME, value=f('LR Pair')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=InteractionMetadataCv.INTERACTION_DIRECTNESS, value=f('Evidence')),
        CV(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=f('Species', map='species_taxon')),
    ),
    membership=MembershipBuilder(
        Member(
            entity=EntityBuilder(
                entity_type=EntityTypeCv.PROTEIN,
                identifiers=IdentifiersBuilder(
                    CV(term=IdentifierNamespaceCv.GENE_NAME_PRIMARY, value=f('Ligand Symbols', extract='primary_gene')),
                    CV(term=IdentifierNamespaceCv.HGNC, value=f('Ligand Species ID', extract='hgnc_id')),
                    CV(term=IdentifierNamespaceCv.ENSEMBL, value=f('Ligand ENSEMBL ID')),
                ),
                annotations=AnnotationsBuilder(
                    CV(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=f('Species', map='species_taxon')),
                    CV(term=lambda row: _location_terms(row.get('Ligand Location'))),
                    CV(term=InterCellAnnotations.LIGAND),
                )
            ),
            annotations=AnnotationsBuilder(
                CV(term=InterCellAnnotations.LIGAND),
            ),
        ),
        Member(
            entity=EntityBuilder(
                entity_type=EntityTypeCv.PROTEIN,
                identifiers=IdentifiersBuilder(
                    CV(term=IdentifierNamespaceCv.GENE_NAME_PRIMARY, value=f('Receptor Symbols', extract='primary_gene')),
                    CV(term=IdentifierNamespaceCv.HGNC, value=f('Receptor Species ID', extract='hgnc_id')),
                    CV(term=IdentifierNamespaceCv.ENSEMBL, value=f('Receptor ENSEMBL ID')),
                ),
                annotations=AnnotationsBuilder(
                    CV(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=f('Species', map='species_taxon')),
                    CV(term=lambda row: _location_terms(row.get('Receptor Location'))),
                    CV(term=InterCellAnnotations.RECEPTOR),
                )
            ),
            annotations=AnnotationsBuilder(
                CV(term=InterCellAnnotations.RECEPTOR),
            ),
        ),
    ),
)

# =============================================================================
# Resource Definition
# =============================================================================

resource = Resource(
    config,
    interactions=Dataset(
        download=download_interactions,
        mapper=interactions_schema,
        raw_parser=iter_csv,
    ),
)
