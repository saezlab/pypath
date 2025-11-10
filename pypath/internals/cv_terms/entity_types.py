"""Entity type controlled vocabularies.

This module contains CV terms for describing entity types (proteins, genes, compounds, etc.).
Identifier namespaces have been moved to identifiers.py.
"""
from .core import CvEnum


class EntityTypeCv(CvEnum):
    """Common PSI-MI and OmniPath entity type terms.

    Defines the types of entities that can be represented in the system,
    including standard PSI-MI types and OmniPath-specific extensions.
    """

    parent_cv_term = "MI:0313"  # interactor type - Molecular species involved in the interaction

    # PSI-MI standard terms
    PROTEIN = "MI:0326"
    GENE = "MI:0250"
    RNA = "MI:0320"
    PROTEIN_COMPLEX = "MI:0315"
    SMALL_MOLECULE = "MI:0328"
    PHENOTYPE = "MI:2261"
    STIMULUS = "MI:2260"

    # OmniPath-specific terms (OM:0010-0099 range)
    PROTEIN_FAMILY = ("OM:0010", "Protein family grouping related proteins")
    LIPID = ("OM:0011", "Lipid molecule or lipid-like compound")
    CV_TERM = ("OM:0012", "Controlled vocabulary term entity")
    INTERACTION = ("OM:0013", "Interaction entity for representing interactions as nodes")
