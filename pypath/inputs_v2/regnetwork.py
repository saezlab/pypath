"""
Parse ConnectomeDB2025 data and emit Entity records.

This module converts ConnectomeDB2025 interactions and complexes into Entity 
records using the declarative schema pattern.
"""

from __future__ import annotations

from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    LicenseCV,
    UpdateCategoryCV,
    ResourceCv,
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
from pypath.inputs_v2.parsers.base import iter_tsv


# =============================================================================
# Resource Configuration
# =============================================================================

config = ResourceConfig(
    id=ResourceCv.REGNETWORK,
    name='RegNetwork',
    url='https://www.zpliulab.cn/RegNetwork/home',
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='40804431',
    primary_category='interactions',
    description=(
        'RegNetwork is an open-source gene regulatory network (GRN) database that'
        'comprehensively curates regulatory relationships among transcription'
        'factors (TFs), microRNAs (miRNAs), and genes in human and mouse.'
        'Establishing a systematic GRN database is pivotal for advancing'
        'researchers’ understanding of gene interactions. Currently, regulatory'
        'interactions are scattered across diverse data resources, posing significant'
        'hurdles for their efficient utilization. Consequently, integrating various'
        'data types and accurately annotating them with confidence scores has become'
        'imperative. In this updated version, we have extensively revised regulatory'
        'relationships among TFs, miRNAs, and genes, introducing novel datasets'
        'involving long noncoding RNAs (lncRNAs) and circular RNAs (circRNAs).'
        'As of now, RegNetwork 2025 comprises 125 319 nodes, including 76 156'
        'for human and 49 163 for mouse, along with 11 107 799 regulatory'
        'interactions, including 7 712 347 for human and 3 395 452 for mouse.'
        'Compared to its predecessor, RegNetwork 2025 has witnessed an over 95%'
        'increase in the total number of regulatory interactions. Furthermore,'
        'we have devised a scoring system to quantify the reliability of regulatory'
        'relationships, enabling us to assemble a more trustworthy core dataset for'
        'users. The incorporation of lncRNA and circRNA has enriched the regulatory'
        'interaction types in RegNetwork 2025 and enhanced the accessibility of'
        'diverse data types. Data are freely available at http://www.zpliulab.cn/RegNetwork/home.'
    ),
)


# =============================================================================
# Download Configurations
# =============================================================================

BASE_URL = 'https://www.zpliulab.cn/RegNetwork/static/data/'

download_human_interactions = Download(
    url=BASE_URL + 'human.7z',
    filename='regnetwork_human_interactions.7z',
    subfolder='regnetwork',
    ext="7z", # TODO: cachedir should support 7z archives
)

download_mouse_interactions = Download(
    url=BASE_URL + 'mouse.7z',
    filename='regnetwork_mouse_interactions.7z',
    subfolder='regnetwork',
    ext="7z", # TODO: cachedir should support 7z archives
)


# =============================================================================
# Processing Helpers
# =============================================================================



# =============================================================================
# Field and Schema Definitions
# =============================================================================

f = FieldConfig()

# -----------------------------------------------------------------------------
# Interactions Schema
# -----------------------------------------------------------------------------

# TODO: add column names to the data files (they were missing)
column_names = ["Regulator Symbol", "Regulator ID", "Target Symbol", "Target ID", "Regulator Type", "Target Type"]

interactions_schema = EntityBuilder(
    entity_type=EntityTypeCv.INTERACTION,
    identifiers=IdentifiersBuilder(),
    annotations=AnnotationsBuilder(),
    membership=MembershipBuilder(
        Member(
            entity=EntityBuilder(
                entity_type=EntityTypeCv.PROTEIN,
                identifiers=IdentifiersBuilder(
                    CV(term=IdentifierNamespaceCv.GENE_NAME_PRIMARY,
                       value=f('Regulator Symbol')),
                    CV(term=IdentifierNamespaceCv.MIRBASE,
                       value=f('Regulator ID')),
                ),
                annotations=AnnotationsBuilder(
                    CV(term=ParticipantMetadataCv.PARTICIPANT_ANNOTATION, value=f('Regulator Type')),
                )
            ),
        ),
        Member(
            entity=EntityBuilder(
                entity_type=EntityTypeCv.GENE,
                identifiers=IdentifiersBuilder(
                    CV(term=IdentifierNamespaceCv.GENE_NAME_PRIMARY,
                       value=f('Target Symbol')),
                    CV(term=IdentifierNamespaceCv.MIRBASE,
                       value=f('Target ID')),
                ),
                annotations=AnnotationsBuilder(
                    CV(term=ParticipantMetadataCv.PARTICIPANT_ANNOTATION, value=f('Target Type')),
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
        raw_parser=iter_tsv,
    ),
    mouse_interactions=Dataset(
        download=download_mouse_interactions,
        mapper=interactions_schema,
        raw_parser=iter_tsv,
    ),
)
