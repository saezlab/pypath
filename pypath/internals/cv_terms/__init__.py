"""Controlled vocabulary terms for the OmniPath build system.

This package provides structured controlled vocabularies (CV) organized by topic,
with support for metadata like definitions and source URLs. All CV terms are based
on either PSI-MI standard accessions or OmniPath-specific OM accessions. PSI-MI
terms only store their accession values here—their explanations and links are
resolved from the MI ontology at runtime.

The CV terms are organized into logical modules:
- entities: Entity types and identifier namespaces
- roles: Biological roles, experimental roles, and membership relationships
- interactions: Interaction types, detection methods, causal effects and mechanisms
- references: Reference types and curation metadata
- licenses: License and update frequency terms

Usage:
    from omnipath_build.utils.cv_terms import EntityTypeCv, BiologicalRoleCv

    # Access as string (accession value)
    print(EntityTypeCv.PROTEIN)  # "MI:0326"

    # Access metadata (available for OmniPath-specific terms)
    print(EntityTypeCv.PROTEIN_FAMILY.definition)  # "Protein family grouping related proteins"
    print(IdentifierNamespaceCv.PDB.source)        # "https://www.rcsb.org"

Accession Ranges:
    PSI-MI standard terms: MI:XXXX
    OmniPath entity types: OM:0010-0099
    OmniPath identifiers: OM:0001-0009 (specialized), OM:0100-0209 (databases/names)
    OmniPath membership roles: OM:0300-0399
    OmniPath update categories: OM:0401-0499
    OmniPath curation: OM:0400-0499
    OmniPath licenses: OM:0501-0599
"""

from typing import Union

from .core import CvEnum
from .entity_types import EntityTypeCv
from .identifiers import IdentifierNamespaceCv
from .annotations import (
    MembershipRoleCv,
    BiologicalRoleCv,
    ExperimentalRoleCv,
    IdentificationMethodCv,
    InteractionTypeCv,
    DetectionMethodCv,
    BiologicalEffectCv,
    CausalMechanismCv,
    PharmacologicalActionCv,
    CausalStatementCv,
    ComplexExpansionCv,
    CurationCv,
    MoleculeAnnotationsCv,
    InteractionParameterCv,
    AffinityUnitCv,
    LigandTypeCv,
    OntologyAnnotationCv,
)
from .resource_metadata import LicenseCV, UpdateCategoryCV
