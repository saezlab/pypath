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
    MoleculeAnnotationsCv,
    MoleculeSubtypeCv,
    UpdateCategoryCV,
    ResourceCv,
    InterCellAnnotations,
)
from pypath.internals.silver_schema import (
    Annotation,
    Entity,
    Identifier,
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
    primary_category='interactions',
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
    ext='csv',
)

download_complexes = Download(
    url=BASE_URL + 'complex_input.csv',
    filename='cellphonedb_complexes.csv',
    subfolder='cellphonedb',
    ext='csv',
)

download_proteins = Download(
    url=BASE_URL + 'protein_input.csv',
    filename='cellphonedb_proteins.csv',
    subfolder='cellphonedb',
    ext='csv',
)


# =============================================================================
# Processing Helpers
# =============================================================================

# Standard UniProt accession regex
UNIPROT_ACC_RE = re.compile(
    r'^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$'
)
HUMAN_TAXON_ID = '9606'
SYNTHETIC_METABOLITE_SYSTEM_RE = re.compile(r'^(.+?)_by[A-Za-z0-9].*')


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


def _synthetic_metabolite_name(val: str) -> str | None:
    """Extract the molecule name from CellPhoneDB synthetic metabolite systems."""
    m = SYNTHETIC_METABOLITE_SYSTEM_RE.match(val or '')
    return m.group(1) if m else None


def _is_synthetic_metabolite_system(val: str) -> bool:
    return _synthetic_metabolite_name(val) is not None


def _extract_partner_name(val: str) -> str | None:
    if UNIPROT_ACC_RE.match(val):
        return None
    return _synthetic_metabolite_name(val) or val


def _extract_synthetic_metabolite_label(val: str) -> str | None:
    return val if _is_synthetic_metabolite_system(val) else None


def _synthetic_metabolite_subtype_term(val: str) -> MoleculeAnnotationsCv | None:
    return (
        MoleculeAnnotationsCv.MOLECULE_SUBTYPE
        if _is_synthetic_metabolite_system(val)
        else None
    )


def _get_partner_type(col: str) -> Any:
    """
    Determine entity type for a partner.

    CellPhoneDB encodes some metabolite systems as names such as
    ``Glutamate_byGLS2_and_SLC1A1``. These are not protein complexes: the
    biological entity is the metabolite, while the suffix is a resource-specific
    label about supporting proteins. We preserve that label as annotation
    instead of turning the proteins into complex members.
    """
    def _type_selector(row: dict[str, Any]) -> EntityTypeCv:
        val = row.get(col, '')
        if _is_synthetic_metabolite_system(val):
            return EntityTypeCv.CHEMICAL
        return (
            EntityTypeCv.PROTEIN 
            if UNIPROT_ACC_RE.match(val) 
            else EntityTypeCv.COMPLEX
        )
    return _type_selector


def _directional_role(row: dict[str, Any], partner: str) -> InterCellAnnotations | None:
    if (row.get('directionality') or '').strip().lower() != 'ligand-receptor':
        return None
    if partner == 'partner_a':
        return InterCellAnnotations.LIGAND
    if partner == 'partner_b':
        return InterCellAnnotations.RECEPTOR
    return None


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
        'partner_name': _extract_partner_name,
        'synthetic_metabolite_label': _extract_synthetic_metabolite_label,
        'synthetic_metabolite_subtype_term': _synthetic_metabolite_subtype_term,
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
        CV(term=IdentifierNamespaceCv.PUBMED, value=f(_source_split, extract='pmid')),
        CV(term=IdentifierNamespaceCv.PUBMED_CENTRAL, value=f(_source_split, extract='pmc')),
    ),
    membership=MembershipBuilder(
        Member(
            entity=EntityBuilder(
                entity_type=_get_partner_type('partner_a'),
                identifiers=IdentifiersBuilder(
                    CV(term=IdentifierNamespaceCv.UNIPROT,
                       value=f('partner_a', extract='uniprot_acc')),
                    CV(term=IdentifierNamespaceCv.NAME,
                       value=f('partner_a', extract='partner_name')),
                ),
                annotations=AnnotationsBuilder(
                    CV(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=HUMAN_TAXON_ID),
                    CV(term=lambda row: _directional_role(row, 'partner_a')),
                    CV(
                        term=f('partner_a', extract='synthetic_metabolite_subtype_term'),
                        value=MoleculeSubtypeCv.METABOLITE,
                    ),
                    CV(
                        term=MoleculeAnnotationsCv.SOURCE_STATUS,
                        value=f('partner_a', extract='synthetic_metabolite_label'),
                    ),
                ),
            ),
            annotations=AnnotationsBuilder(
                CV(term=lambda row: _directional_role(row, 'partner_a')),
            ),
        ),
        Member(
            entity=EntityBuilder(
                entity_type=_get_partner_type('partner_b'),
                identifiers=IdentifiersBuilder(
                    CV(term=IdentifierNamespaceCv.UNIPROT,
                       value=f('partner_b', extract='uniprot_acc')),
                    CV(term=IdentifierNamespaceCv.NAME, 
                       value=f('partner_b', extract='partner_name')),
                ),
                annotations=AnnotationsBuilder(
                    CV(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=HUMAN_TAXON_ID),
                    CV(term=lambda row: _directional_role(row, 'partner_b')),
                    CV(
                        term=f('partner_b', extract='synthetic_metabolite_subtype_term'),
                        value=MoleculeSubtypeCv.METABOLITE,
                    ),
                    CV(
                        term=MoleculeAnnotationsCv.SOURCE_STATUS,
                        value=f('partner_b', extract='synthetic_metabolite_label'),
                    ),
                ),
            ),
            annotations=AnnotationsBuilder(
                CV(term=lambda row: _directional_role(row, 'partner_b')),
            ),
        ),
    ),
)

# -----------------------------------------------------------------------------
# Complexes Schema
# -----------------------------------------------------------------------------

protein_complexes_schema = EntityBuilder(
    entity_type=EntityTypeCv.COMPLEX,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.NAME, value=f('complex_name')),
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
            entity_annotations=AnnotationsBuilder(
                CV(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=HUMAN_TAXON_ID),
            ),
        )
    ),
)


def _identifier(type_: object, value: object) -> Identifier | None:
    value = str(value or '').strip()
    return Identifier(type=type_, value=value) if value else None


def _annotation(term: object, value: object = None) -> Annotation | None:
    value = str(value or '').strip() if value is not None else None
    return Annotation(term=term, value=value) if value else None


def _identifiers(*items: Identifier | None) -> list[Identifier]:
    out: list[Identifier] = []
    seen: set[tuple[object, str]] = set()
    for item in items:
        if item is None:
            continue
        key = (item.type, item.value)
        if key in seen:
            continue
        out.append(item)
        seen.add(key)
    return out


def _annotations(*items: Annotation | None) -> list[Annotation] | None:
    out: list[Annotation] = []
    seen: set[tuple[object, object, object]] = set()
    for item in items:
        if item is None:
            continue
        key = (item.term, item.value, item.units)
        if key in seen:
            continue
        out.append(item)
        seen.add(key)
    return out or None


def complexes_schema(row: dict[str, Any]) -> Entity:
    name = str(row.get('complex_name') or '').strip()
    metabolite_name = _synthetic_metabolite_name(name)
    if not metabolite_name:
        return protein_complexes_schema(row)

    return Entity(
        type=EntityTypeCv.CHEMICAL,
        identifiers=_identifiers(
            _identifier(IdentifierNamespaceCv.NAME, metabolite_name),
        ),
        annotations=_annotations(
            Annotation(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=HUMAN_TAXON_ID),
            Annotation(
                term=MoleculeAnnotationsCv.MOLECULE_SUBTYPE,
                value=MoleculeSubtypeCv.METABOLITE,
            ),
            _annotation(MoleculeAnnotationsCv.SOURCE_STATUS, name),
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
