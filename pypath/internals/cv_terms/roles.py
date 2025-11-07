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

    # Ontological and membership relationships (OM:0300-0399 range)
    IS_A = ("OM:0300", "Hierarchical parent-child relationship in ontologies")
    PART_OF = ("OM:0301", "Component or part relationship")
    MEMBER_OF = ("OM:0302", "Generic membership in a collection or group")
    IS_ANNOTATED_AS = ("OM:0303", "Annotation relationship linking entity to CV term")


class BiologicalRoleCv(CvEnum):
    """Biological role terms for interaction participants.

    Describes the functional role of an entity in a biological process
    or interaction (e.g., enzyme, substrate, inhibitor).
    """

    # PSI-MI standard biological roles
    ENZYME = ("MI:0501", "Catalyst of a biochemical reaction")
    SUBSTRATE = ("MI:0502", "Molecule acted upon by an enzyme")
    INHIBITOR = ("MI:0586", "Molecule that decreases activity of a target")
    STIMULATOR = ("MI:0840", "Molecule that increases activity of a target")
    ALLOSTERIC_EFFECTOR = ("MI:1160", "Molecule that binds allosterically to modulate activity")
    REGULATOR_TARGET = ("MI:2275", "Target of a regulatory interaction")


class ExperimentalRoleCv(CvEnum):
    """Experimental role terms for interaction participants.

    Describes the experimental context or role in which an entity
    was identified in a study (e.g., bait, prey, control).
    """

    # PSI-MI standard experimental roles
    BAIT = ("MI:0496", "Protein used as bait in pull-down or affinity experiments")
    PREY = ("MI:0498", "Protein identified as binding partner in pull-down experiments")
    NEUTRAL_COMPONENT = ("MI:0497", "Component with no specific experimental bias")
    UNSPECIFIED_ROLE = ("MI:0499", "Role not specified or not applicable")


class IdentificationMethodCv(CvEnum):
    """Participant identification method terms.

    Describes the method used to identify or detect a participant
    in an interaction study.
    """

    # PSI-MI standard identification methods
    GENERIC = ("MI:0002", "Generic or unspecified identification method")
    MASS_SPECTROMETRY = ("MI:0427", "Mass spectrometry-based identification")
    SEQUENCE_TAG = ("MI:0102", "Identification by sequence tag or epitope tag")
