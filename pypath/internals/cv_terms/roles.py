"""Role-related controlled vocabularies.

This module contains CV terms for describing roles of participants in interactions
and memberships, including biological roles, experimental roles, identification
methods, and membership relationships.
"""
from .core import CvEnum


class MembershipRoleCv(CvEnum):
    """Relationship types for entity membership (is_member_of field).

    Defines how entities relate to their parent entities in hierarchical
    or compositional relationships (e.g., proteins in complexes, ontology terms).
    """

    parent_cv_term = ("OM:0300", "Membership role term - Describes the membership relationship of an entity.")  # OmniPath-specific term - no PSI-MI parent

    # Ontological and membership relationships (OM:0300-0399 range)
    IS_A = ("OM:0301", "Hierarchical parent-child relationship in ontologies")
    PART_OF = ("OM:0302", "Component or part relationship")
    MEMBER_OF = ("OM:0303", "Generic membership in a collection or group")
    IS_ANNOTATED_AS = ("OM:0304", "Annotation relationship linking entity to CV term")


class BiologicalRoleCv(CvEnum):
    """Biological role terms for interaction participants.

    Describes the functional role of an entity in a biological process
    or interaction (e.g., enzyme, substrate, inhibitor).
    """

    parent_cv_term = "MI:0500"  # biological role - Physiological role of an interactor in a cell or in vivo environment

    # PSI-MI standard biological roles
    ENZYME = "MI:0501"
    SUBSTRATE = "MI:0502"
    INHIBITOR = "MI:0586"
    STIMULATOR = "MI:0840"
    ALLOSTERIC_EFFECTOR = "MI:1160"
    REGULATOR_TARGET = "MI:2275"


class ExperimentalRoleCv(CvEnum):
    """Experimental role terms for interaction participants.

    Describes the experimental context or role in which an entity
    was identified in a study (e.g., bait, prey, control).
    """

    parent_cv_term = "MI:0495"  # experimental role - Role played by the participant within the experiment

    # PSI-MI standard experimental roles
    BAIT = "MI:0496"
    PREY = "MI:0498"
    NEUTRAL_COMPONENT = "MI:0497"
    UNSPECIFIED_ROLE = "MI:0499"


class IdentificationMethodCv(CvEnum):
    """Participant identification method terms.

    Describes the method used to identify or detect a participant
    in an interaction study.
    """

    parent_cv_term = "MI:0002"  # participant identification method - Method to determine the molecules involved in the interaction

    # PSI-MI standard identification methods
    GENERIC = "MI:0002"
    MASS_SPECTROMETRY = "MI:0427"
    SEQUENCE_TAG = "MI:0102"
