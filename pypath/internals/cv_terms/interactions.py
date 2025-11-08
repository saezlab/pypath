"""Interaction-related controlled vocabularies.

This module contains CV terms for describing interactions, including interaction types,
detection methods, biological effects, causal mechanisms, causal statements, and
complex expansion methods.
"""
from .core import CvEnum


class InteractionTypeCv(CvEnum):
    """Interaction type terms from PSI-MI.

    Describes the nature or type of relationship between interacting entities
    (e.g., physical binding, functional association, enzymatic reaction).
    """

    parent_cv_term = "MI:0190"  # interaction type - Connection between molecule

    # PSI-MI standard interaction types
    COLOCALIZATION = "MI:0403"
    FUNCTIONAL_ASSOCIATION = "MI:2286"
    PHYSICAL_ASSOCIATION = "MI:0915"
    DIRECT_INTERACTION = "MI:0407"
    PHOSPHORYLATION_REACTION = "MI:0217"
    PHENOTYPE_RESULT = "MI:2283"


class DetectionMethodCv(CvEnum):
    """Experimental detection method terms from PSI-MI.

    Describes the experimental technique used to detect or validate
    an interaction between entities.
    """

    parent_cv_term = "MI:0001"  # interaction detection method - Method to determine the interaction

    # PSI-MI standard detection methods
    AFFINITY_CHROMATOGRAPHY = "MI:0004"
    COIMMUNOPRECIPITATION = "MI:0019"
    PULL_DOWN = "MI:0096"
    INFERRED_BY_CURATOR = "MI:0364"


class BiologicalEffectCv(CvEnum):
    """Biological effect terms describing causal outcomes.

    Describes the functional consequence of an interaction on the
    target entity (e.g., increase or decrease in activity or quantity).
    """

    parent_cv_term = "MI:2233"  # causal interaction - Binary causative relationships between biological entities

    # PSI-MI standard biological effects
    UP_REGULATES_ACTIVITY = "MI:2236"
    DOWN_REGULATES_ACTIVITY = "MI:2241"
    UP_REGULATES_QUANTITY = "MI:2237"
    DOWN_REGULATES_QUANTITY = "MI:2242"


class CausalMechanismCv(CvEnum):
    """Causal mechanism terms from PSI-MI.

    Describes the molecular mechanism by which a causal effect is achieved
    (e.g., transcriptional regulation, post-translational modification).
    """

    parent_cv_term = "MI:2233"  # causal interaction - Binary causative relationships between biological entities

    # PSI-MI standard causal mechanisms
    TRANSCRIPTIONAL_REGULATION = "MI:2247"
    TRANSLATION_REGULATION = "MI:2248"
    POST_TRANSLATIONAL_REGULATION = "MI:2249"


class CausalStatementCv(CvEnum):
    """Causal statement terms from PSI-MI.

    Provides detailed causal statements combining direction (up/down regulation)
    with the specific mechanism or target (activity, quantity, stability, expression).
    """

    parent_cv_term = "MI:2234"  # causal statement - The effect of modulator entity A on a modulated entity B

    # Down-regulation terms
    DOWN_REGULATES = "MI:2240"
    DOWN_REGULATES_ACTIVITY = "MI:2241"
    DOWN_REGULATES_QUANTITY = "MI:2242"
    DOWN_REGULATES_QUANTITY_BY_DESTABLIZATION = "MI:2244"
    DOWN_REGULATES_QUANTITY_BY_REPRESSION = "MI:2243"

    # Up-regulation terms
    UP_REGULATES = "MI:2235"
    UP_REGULATES_ACTIVITY = "MI:2236"
    UP_REGULATES_QUANTITY = "MI:2237"
    UP_REGULATES_QUANTITY_BY_EXPRESSION = "MI:2238"
    UP_REGULATES_QUANTITY_BY_STABILIZATION = "MI:2239"


class ComplexExpansionCv(CvEnum):
    """Complex expansion strategy terms from PSI-MI.

    Describes the method used to expand n-ary complexes into binary interactions
    for analysis and representation purposes.
    """

    parent_cv_term = "MI:1059"  # complex expansion - The method by which complex n-ary data is expanded into binary data

    # PSI-MI standard complex expansion methods
    BIPARTITE_EXPANSION = "MI:1062"
    MATRIX_EXPANSION = "MI:1061"
    SPOKE_EXPANSION = "MI:1060"
