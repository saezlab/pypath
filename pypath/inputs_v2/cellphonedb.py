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

# Standard UniProt accession regex
UNIPROT_ACC_RE = re.compile(
    r'^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$'
)


def _extract_pmid(token: str) -> str | None:
    """Extract PubMed ID from a token."""
    m = re.search(r'PMID:?\s*(\d+)', token, re.IGNORECASE)
    return m.group(1) if m else None


def _extract_pmc(token: str) -> str | None:
    """Extract PubMed Central ID from a token."""
    m = re.search(r'PMC\s*(\d+)', token, re.IGNORECASE)
    return f'PMC{m.group(1)}' if m else None


def _extract_comment(token: str) -> str | None:
    """Return the token if it's not a PMID or PMC."""
    if re.search(r'PMID|PMC', token, re.IGNORECASE):
        return None
    return token.strip()


def _source_split(row: dict[str, Any]) -> list[str]:
    """Split the source column into individual tokens using regex."""
    return re.split(r'[;,]\s*', row.get('source') or '')


def _extract_uniprot_acc(val: str) -> str | None:
    """Return the value if it matches the UniProt accession pattern."""
    return val if UNIPROT_ACC_RE.match(val) else None


def _extract_non_uniprot(val: str) -> str | None:
    """Return the value if it does NOT match the UniProt accession pattern."""
    return val if not UNIPROT_ACC_RE.match(val) else None


def _get_partner_type(col: str) -> Any:
    """
    Determine entity type for a partner (Protein or Complex).
    """
    def _type_selector(row: dict[str, Any]) -> EntityTypeCv:
        val = row.get(col, '')
        return (
            EntityTypeCv.PROTEIN 
            if UNIPROT_ACC_RE.match(val) 
            else EntityTypeCv.COMPLEX
        )
    return _type_selector


# =============================================================================
# Field and Schema Definitions
# =============================================================================

f = FieldConfig(
    extract={
        'pmid': _extract_pmid,
        'pmc': _extract_pmc,
        'comment': _extract_comment,
        'uniprot_acc': _extract_uniprot_acc,
        'non_uniprot': _extract_non_uniprot,
    },
)

# -----------------------------------------------------------------------------
# Interactions Schema
# -----------------------------------------------------------------------------

interactions_schema = EntityBuilder(
    entity_type=EntityTypeCv.INTERACTION,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.NAME, value=f('interactors')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=InteractionMetadataCv.INTERACTION_ANNOTATION, value=f('classification')),
        CV(term=InteractionMetadataCv.INTERACTION_ANNOTATION, value=f('directionality')),
        CV(term=IdentifierNamespaceCv.PUBMED, value=f(_source_split, extract='pmid')),
        CV(term=IdentifierNamespaceCv.PUBMED_CENTRAL, value=f(_source_split, extract='pmc')),
        CV(term=CurationCv.COMMENT, value=f(_source_split, extract='comment')),
        CV(term=CurationCv.COMMENT, value=f('version')),
    ),
    membership=MembershipBuilder(
        Member(
            entity=EntityBuilder(
                entity_type=_get_partner_type('partner_a'),
                identifiers=IdentifiersBuilder(
                    CV(term=IdentifierNamespaceCv.UNIPROT, value=f('partner_a', extract='uniprot_acc')),
                    CV(term=IdentifierNamespaceCv.NAME, value=f('partner_a', extract='non_uniprot')),
                ),
            ),
        ),
        Member(
            entity=EntityBuilder(
                entity_type=_get_partner_type('partner_b'),
                identifiers=IdentifiersBuilder(
                    CV(term=IdentifierNamespaceCv.UNIPROT, value=f('partner_b', extract='uniprot_acc')),
                    CV(term=IdentifierNamespaceCv.NAME, value=f('partner_b', extract='non_uniprot')),
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
                    value=f(
                        lambda row: [
                            row.get(f'uniprot_{i}')
                            for i in range(1, 5)
                            if row.get(f'uniprot_{i}')
                        ],
                    ),
                ),
            ),
        )
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
)
