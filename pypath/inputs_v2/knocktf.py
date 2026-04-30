"""
Parse CellPhoneDB data and emit Entity records.

This module converts CellPhoneDB interactions and complexes into Entity 
records using the declarative schema pattern.
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
from pypath.inputs_v2.parsers.base import iter_tsv


# =============================================================================
# Resource Configuration
# =============================================================================

config = ResourceConfig(
    id=ResourceCv.KNOCKTF,
    name='KnockTF 2.0',
    url='http://www.licpathway.net:8081/KnockTFv2/',
    license=LicenseCV.CC_BY_NC_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='37956336',
    primary_category='interactions',
    description=(
        'KnockTF 2.0 aims to provide comprehensive T(co)F knockdown/knockout'
        'dataset resource across multiple tissue/cell types of different'
        'species. The current version of KnockTF stores 1468 manually'
        'curated RNA-seq and microarray datasets associated with 612 TFs'
        'and 172 TcoFs disrupted by different knockdown/knockout techniques'
        'and across multiple tissue/cell types in humans, mice, Arabidopsis'
        'thaliana and Zea mays.'
        '\n'
        'KnockTF 2.0 not only provides comprehensive gene expression information for'
        'T(co)F target genes of interest, but also collects upstream pathway information'
        'for T(co)Fs and various functional annotation information for downstream target'
        'genes, such as GO/KEGG pathway enrichment, hierarchical clustering analysis and'
        'differentially expressed analysis. Furthermore, KnockTF 2.0 provides the detailed'
        'and abundant (epi)genetic annotation information for T(co)F target genes, including'
        'super-enhancers, enhancers, transcription factor binding sites (TFBSs), common SNPs,'
        'risk SNPs, LD SNPs, expression quantitative trait locus (eQTL), methylation sites,'
        'DNase I hypersensitivity sites (DHSs), chromatin interactions, chromatin accessibility'
        'regions, CRISPR/Cas9 target sites, and topologically associating domains (TADs).'
        '\n'
        'KnockTF 2.0 will facilitate the identification of functional T(co)Fs and target genes,'
        'and benefit the investigation of their roles in the physiological and pathological processes.'
    ),
)


# =============================================================================
# Download Configurations
# =============================================================================

BASE_URL = 'http://www.licpathway.net:8081/KnockTFv2/public/download_anno/'

download_DE_human = Download(
    url=BASE_URL + 'knocktf_v2_main_human.txt',
    filename='knocktf_v2_main_human.txt',
    subfolder='knocktf',
)

download_DE_mouse = Download(
    url=BASE_URL + 'knocktf_v2_main_mouse.txt',
    filename='knocktf_v2_main_mouse.txt',
    subfolder='knocktf',
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

interactions_schema = EntityBuilder(
    entity_type=EntityTypeCv.INTERACTION,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.KNOCKTF, value=f('Sample_ID')), # not an identifier for the row
    ),
    annotations=AnnotationsBuilder( # TODO: specify term
        CV(term=InteractionMetadataCv.INTERACTION_ANNOTATION, value=f('Mean_Case')),
        CV(term=InteractionMetadataCv.INTERACTION_ANNOTATION, value=f('Mean_Control')),
        CV(term=InteractionMetadataCv.INTERACTION_ANNOTATION, value=f('FC')),
        CV(term=InteractionMetadataCv.INTERACTION_ANNOTATION, value=f('Log2FC')),
        CV(term=InteractionMetadataCv.INTERACTION_ANNOTATION, value=f('Rank')),
        CV(term=InteractionMetadataCv.INTERACTION_ANNOTATION, value=f('P_value')),
        CV(term=InteractionMetadataCv.INTERACTION_ANNOTATION, value=f('Corrected_P')),
        CV(term=InteractionMetadataCv.INTERACTION_ANNOTATION, value=f('up_down')),
        CV(term=InteractionMetadataCv.INTERACTION_ANNOTATION, value=f('TE_TF')),
        CV(term=InteractionMetadataCv.INTERACTION_ANNOTATION, value=f('promoter_TF')),
        CV(term=InteractionMetadataCv.INTERACTION_ANNOTATION, value=f('SE_TF')),
        CV(term=InteractionMetadataCv.INTERACTION_ANNOTATION, value=f('TE_TF_name')),
        CV(term=InteractionMetadataCv.INTERACTION_ANNOTATION, value=f('promoter_TF_name')),
        CV(term=InteractionMetadataCv.INTERACTION_ANNOTATION, value=f('SE_TF_name')),
    ),
    membership=MembershipBuilder(
        Member(
            entity=EntityBuilder(
                entity_type=EntityTypeCv.PROTEIN,
                identifiers=IdentifiersBuilder(
                    CV(term=IdentifierNamespaceCv.GENE_NAME_PRIMARY,
                       value=f('TF')),
                ),
            ),
        ),
        Member(
            entity=EntityBuilder(
                entity_type=EntityTypeCv.GENE,
                identifiers=IdentifiersBuilder(
                    CV(term=IdentifierNamespaceCv.GENE_NAME_PRIMARY,
                       value=f('NCBI_Official_Symbol')),
                    CV(term=IdentifierNamespaceCv.GENE_NAME_SYNONYM,
                       value=f('Gene')),
                    CV(term=IdentifierNamespaceCv.ENTREZ, 
                       value=f('Entrez_ID')),
                    CV(term=IdentifierNamespaceCv.ENSEMBL, 
                       value=f('Ensembl_ID')),
                ),
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
        download=download_DE_human,
        mapper=interactions_schema,
        raw_parser=iter_tsv,
    ),
    mouse_interactions=Dataset(
        download=download_DE_mouse,
        mapper=interactions_schema,
        raw_parser=iter_tsv,
    ),
)
