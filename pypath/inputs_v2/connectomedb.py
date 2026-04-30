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
    LicenseCV,
    UpdateCategoryCV,
    ResourceCv,
    CurationCv,
    InteractionMetadataCv,
    ParticipantMetadataCv,
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
    update_category=UpdateCategoryCV.REGULAR,
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

download_human_interactions = Download(
    url=BASE_URL + 'ConnectomeDB2025_human.csv',
    filename='connectomedb_human_interactions.csv',
    subfolder='connectomedb2025',
)

download_mouse_interactions = Download(
    url=BASE_URL + 'ConnectomeDB2025_mouse.csv',
    filename='connectomedb_mouse_interactions.csv',
    subfolder='connectomedb2025',
)


# =============================================================================
# Processing Helpers
# =============================================================================

_symbols_pat = re.compile(r"^(\w+)(?:\s*\((.+)\))?")

def _extract_primary_gene(token: str):
    return _symbols_pat.search(token).group(1)

def _extract_gene_alias(token: str):
    return _symbols_pat.search(token).group(2).replace(", ", ";")

_cdb_pat = re.compile(r"^CDB\d{2}:(\d+)", re.IGNORECASE)

def _extract_cdb(val: str) -> str | None:
    return val if _cdb_pat.match(val) else None



# =============================================================================
# Field and Schema Definitions
# =============================================================================

f = FieldConfig(
    extract={
        'primary_gene': _extract_primary_gene,
        'gene_alias': _extract_gene_alias,
        'cdb': _extract_cdb,
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
        CV(term=InteractionMetadataCv.INTERACTION_ANNOTATION, value=f('Evidence')),
        CV(term=CurationCv.COMMENT, value=f('AI summary')),
    ),
    membership=MembershipBuilder(
        Member(
            entity=EntityBuilder(
                entity_type=EntityTypeCv.PROTEIN,
                identifiers=IdentifiersBuilder(
                    CV(term=IdentifierNamespaceCv.GENE_NAME_PRIMARY,
                       value=f('Ligand Symbols', extract='primary_gene')),
                ),
                annotations=AnnotationsBuilder(
                    CV(term=IdentifierNamespaceCv.HGNC, value=f('Ligand HGNC ID')),
                    CV(term=IdentifierNamespaceCv.ENSEMBL, value=f('Ligand ENSEMBL ID')),
                    CV(term=ParticipantMetadataCv.ALIAS, value=f('Ligand Symbols', extract='gene_alias')),
                    CV(term=ParticipantMetadataCv.PARTICIPANT_ANNOTATION, value=f('Ligand Location')),
                )
            ),
        ),
        Member(
            entity=EntityBuilder(
                entity_type=EntityTypeCv.PROTEIN,
                identifiers=IdentifiersBuilder(
                    CV(term=IdentifierNamespaceCv.GENE_NAME_PRIMARY,
                       value=f('Receptor Symbols', extract='primary_gene')),
                ),
                annotations=AnnotationsBuilder(
                    CV(term=IdentifierNamespaceCv.HGNC, value=f('Receptor HGNC ID')),
                    CV(term=IdentifierNamespaceCv.ENSEMBL, value=f('Receptor ENSEMBL ID')),
                    CV(term=ParticipantMetadataCv.ALIAS, value=f('Receptor Symbols', extract='gene_alias')),
                    CV(term=ParticipantMetadataCv.PARTICIPANT_ANNOTATION, value=f('Receptor Location')),
                )
            ),
        ),
    ),
)

# =============================================================================
# Resource Definition
# =============================================================================

resource = Resource(
    config,
    human_interactions=Dataset(
        download=download_human_interactions,
        mapper=interactions_schema,
        raw_parser=iter_csv,
    ),
    mouse_interactions=Dataset(
        download=download_mouse_interactions,
        mapper=interactions_schema,
        raw_parser=iter_csv,
    ),
)
