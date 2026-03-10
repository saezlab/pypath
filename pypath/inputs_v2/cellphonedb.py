"""
Parse CellPhoneDB data and emit Entity records.

This module converts CellPhoneDB interactions, complexes, and protein 
annotations into Entity records using the declarative schema pattern.
"""

from __future__ import annotations

import re
from typing import Any

from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    LicenseCV,
    UpdateCategoryCV,
    ResourceCv,
    CurationCv,
    InteractionMetadataCv,
)
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    Member,
    MembershipBuilder,
    MembersFromList,
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.base import iter_csv


# =============================================================================
# Resource Configuration
# =============================================================================

config = ResourceConfig(
    id=ResourceCv.CELLPHONEDB,
    name='CellPhoneDB',
    url='https://www.cellphonedb.org/',
    license=LicenseCV.MIT,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='40133495',
    description=(
        'CellPhoneDB is a publicly available repository of curated receptors, '
        'ligands and their interactions, designed to enable the analysis of '
        'cell-cell communication from single-cell transcriptomics data.'
    ),
)


# =============================================================================
# Download Configurations
# =============================================================================

BASE_URL = 'https://raw.githubusercontent.com/ventolab/cellphonedb-data/master/data/'

download_interactions = Download(
    url=BASE_URL + 'interaction_input.csv',
    filename='cellphonedb_interactions.csv',
    subfolder='cellphonedb',
)

download_complexes = Download(
    url=BASE_URL + 'complex_input.csv',
    filename='cellphonedb_complexes.csv',
    subfolder='cellphonedb',
)

download_proteins = Download(
    url=BASE_URL + 'protein_input.csv',
    filename='cellphonedb_proteins.csv',
    subfolder='cellphonedb',
)


# =============================================================================
# Processing Helpers
# =============================================================================

SOURCE_MAP = {
    'uniprot': IdentifierNamespaceCv.UNIPROT,
    'reactome': IdentifierNamespaceCv.REACTOME_ID,
    'iuphar': IdentifierNamespaceCv.GUIDETOPHARMA,
    'guidetopharmacology.org': IdentifierNamespaceCv.GUIDETOPHARMA,
    'intact': IdentifierNamespaceCv.INTACT,
    'mint': IdentifierNamespaceCv.MINT,
    'dip': IdentifierNamespaceCv.DIP,
    'bind': IdentifierNamespaceCv.BIND,
    'imex': IdentifierNamespaceCv.IMEX,
}


def extract_sources(row: dict[str, Any]) -> list[tuple[IdentifierNamespaceCv, str]]:
    """
    Extract multiple sources from the CellPhoneDB source column.
    """
    source_str = row.get('source', '')
    if not source_str:
        return []

    results = []
    
    # Split by common delimiters
    tokens = re.split(r'[;,]\s*', source_str)
    
    for token in tokens:
        token = token.strip()
        if not token:
            continue
            
        # Handle PMIDs
        pmid_match = re.search(r'PMID:?\s*(\d+)', token, re.IGNORECASE)
        if pmid_match:
            results.append((IdentifierNamespaceCv.PUBMED, pmid_match.group(1)))
            continue
            
        # Handle PMCs
        pmc_match = re.search(r'PMC\s*(\d+)', token, re.IGNORECASE)
        if pmc_match:
            results.append((IdentifierNamespaceCv.PUBMED_CENTRAL, f'PMC{pmc_match.group(1)}'))
            continue

        # Handle known string sources
        token_lower = token.lower()
        if token_lower in SOURCE_MAP:
            results.append((SOURCE_MAP[token_lower], token))
            
    return results


# =============================================================================
# Field and Schema Definitions
# =============================================================================

f = FieldConfig()

# -----------------------------------------------------------------------------
# Interactions Schema
# -----------------------------------------------------------------------------

interactions_schema = EntityBuilder(
    entity_type=EntityTypeCv.INTERACTION,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.CELLPHONEDB, value=f('id_cp_interaction')),
        CV(term=IdentifierNamespaceCv.NAME, value=f('interactors')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=InteractionMetadataCv.INTERACTION_ANNOTATION, value=f('classification')),
        CV(term=InteractionMetadataCv.CONTROL_TYPE, value=f('directionality')),
        CV(term=IdentifierNamespaceCv.PUBMED, 
           value=lambda row: [v for t, v in extract_sources(row) if t == IdentifierNamespaceCv.PUBMED]),
        CV(term=IdentifierNamespaceCv.PUBMED_CENTRAL,
            value=lambda row: [v for t, v in extract_sources(row) if t == IdentifierNamespaceCv.PUBMED_CENTRAL]),
        CV(term=IdentifierNamespaceCv.UNIPROT, 
           value=lambda row: [v for t, v in extract_sources(row) if t == IdentifierNamespaceCv.UNIPROT]),
        CV(term=IdentifierNamespaceCv.REACTOME_ID, 
           value=lambda row: [v for t, v in extract_sources(row) if t == IdentifierNamespaceCv.REACTOME_ID]),
        CV(term=IdentifierNamespaceCv.GUIDETOPHARMA, 
           value=lambda row: [v for t, v in extract_sources(row) if t == IdentifierNamespaceCv.GUIDETOPHARMA]),
        CV(term=IdentifierNamespaceCv.INTACT, 
           value=lambda row: [v for t, v in extract_sources(row) if t == IdentifierNamespaceCv.INTACT]),
        CV(term=IdentifierNamespaceCv.MINT, 
           value=lambda row: [v for t, v in extract_sources(row) if t == IdentifierNamespaceCv.MINT]),
        CV(term=IdentifierNamespaceCv.DIP, 
           value=lambda row: [v for t, v in extract_sources(row) if t == IdentifierNamespaceCv.DIP]),
        CV(term=IdentifierNamespaceCv.BIND, 
           value=lambda row: [v for t, v in extract_sources(row) if t == IdentifierNamespaceCv.BIND]),
        CV(term=IdentifierNamespaceCv.IMEX, 
           value=lambda row: [v for t, v in extract_sources(row) if t == IdentifierNamespaceCv.IMEX]),
        CV(term=CurationCv.COMMENT, value=f('version')),
    ),
    membership=MembershipBuilder(
        Member(
            entity=EntityBuilder(
                entity_type=EntityTypeCv.PROTEIN, 
                identifiers=IdentifiersBuilder(
                    CV(term=IdentifierNamespaceCv.UNIPROT, value=f('partner_a')),
                    CV(term=IdentifierNamespaceCv.NAME, value=f('partner_a')),
                ),
            ),
        ),
        Member(
            entity=EntityBuilder(
                entity_type=EntityTypeCv.PROTEIN,
                identifiers=IdentifiersBuilder(
                    CV(term=IdentifierNamespaceCv.UNIPROT, value=f('partner_b')),
                    CV(term=IdentifierNamespaceCv.NAME, value=f('partner_b')),
                ),
            ),
        ),
    ),
)

# -----------------------------------------------------------------------------
# Complexes Schema
# -----------------------------------------------------------------------------

complexes_schema = EntityBuilder(
    entity_type=EntityTypeCv.COMPLEX,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.NAME, value=f('complex_name')),
        CV(term=IdentifierNamespaceCv.CELLPHONEDB, value=f('complex_name')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=CurationCv.COMMENT, value=f('version')),
    ),
    membership=MembershipBuilder(
        MembersFromList(
            entity_type=EntityTypeCv.PROTEIN,
            identifiers=IdentifiersBuilder(
                CV(
                    term=IdentifierNamespaceCv.UNIPROT,
                    value=lambda row: [row.get(f'uniprot_{i}') for i in range(1, 5) if row.get(f'uniprot_{i}')],
                ),
            ),
        )
    ),
)

# -----------------------------------------------------------------------------
# Proteins Schema
# -----------------------------------------------------------------------------

proteins_schema = EntityBuilder(
    entity_type=EntityTypeCv.PROTEIN,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.UNIPROT, value=f('uniprot')),
        CV(term=IdentifierNamespaceCv.NAME, value=f('protein_name')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=CurationCv.COMMENT, value=f('version')),
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
    complexes=Dataset(
        download=download_complexes,
        mapper=complexes_schema,
        raw_parser=iter_csv,
    ),
    proteins=Dataset(
        download=download_proteins,
        mapper=proteins_schema,
        raw_parser=iter_csv,
    ),
)
