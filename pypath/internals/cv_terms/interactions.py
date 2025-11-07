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

    # PSI-MI standard interaction types
    COLOCALIZATION = ("MI:0403", "Entities found in the same cellular location")
    FUNCTIONAL_ASSOCIATION = ("MI:2286", "Functional relationship without direct physical contact")
    PHYSICAL_ASSOCIATION = ("MI:0915", "Physical interaction or binding")
    DIRECT_INTERACTION = ("MI:0407", "Direct physical interaction between entities")
    PHOSPHORYLATION_REACTION = ("MI:0217", "Phosphorylation enzymatic reaction")
    PHENOTYPE_RESULT = ("MI:2283", "Interaction resulting in a phenotypic effect")


class DetectionMethodCv(CvEnum):
    """Experimental detection method terms from PSI-MI.

    Describes the experimental technique used to detect or validate
    an interaction between entities.
    """

    # PSI-MI standard detection methods
    AFFINITY_CHROMATOGRAPHY = ("MI:0004", "Affinity chromatography technique")
    COIMMUNOPRECIPITATION = ("MI:0019", "Co-immunoprecipitation experiment")
    PULL_DOWN = ("MI:0096", "Pull-down assay")
    INFERRED_BY_CURATOR = ("MI:0364", "Interaction inferred by manual curation")


class BiologicalEffectCv(CvEnum):
    """Biological effect terms describing causal outcomes.

    Describes the functional consequence of an interaction on the
    target entity (e.g., increase or decrease in activity or quantity).
    """

    # PSI-MI standard biological effects
    UP_REGULATES_ACTIVITY = ("MI:2236", "Increases the activity of the target")
    DOWN_REGULATES_ACTIVITY = ("MI:2241", "Decreases the activity of the target")
    UP_REGULATES_QUANTITY = ("MI:2237", "Increases the quantity or abundance of the target")
    DOWN_REGULATES_QUANTITY = ("MI:2242", "Decreases the quantity or abundance of the target")


class CausalMechanismCv(CvEnum):
    """Causal mechanism terms from PSI-MI.

    Describes the molecular mechanism by which a causal effect is achieved
    (e.g., transcriptional regulation, post-translational modification).
    """

    # PSI-MI standard causal mechanisms
    TRANSCRIPTIONAL_REGULATION = ("MI:2247", "Regulation at the transcriptional level")
    TRANSLATION_REGULATION = ("MI:2248", "Regulation at the translational level")
    POST_TRANSLATIONAL_REGULATION = ("MI:2249", "Regulation by post-translational modification")


class CausalStatementCv(CvEnum):
    """Causal statement terms from PSI-MI.

    Provides detailed causal statements combining direction (up/down regulation)
    with the specific mechanism or target (activity, quantity, stability, expression).
    """

    # Down-regulation terms
    DOWN_REGULATES = ("MI:2240", "Decreases or inhibits the target (generic)")
    DOWN_REGULATES_ACTIVITY = ("MI:2241", "Decreases the activity of the target")
    DOWN_REGULATES_QUANTITY = ("MI:2242", "Decreases the quantity of the target")
    DOWN_REGULATES_QUANTITY_BY_DESTABLIZATION = (
        "MI:2244",
        "Decreases quantity by destabilizing the target"
    )
    DOWN_REGULATES_QUANTITY_BY_REPRESSION = (
        "MI:2243",
        "Decreases quantity by repressing expression"
    )

    # Up-regulation terms
    UP_REGULATES = ("MI:2235", "Increases or activates the target (generic)")
    UP_REGULATES_ACTIVITY = ("MI:2236", "Increases the activity of the target")
    UP_REGULATES_QUANTITY = ("MI:2237", "Increases the quantity of the target")
    UP_REGULATES_QUANTITY_BY_EXPRESSION = (
        "MI:2238",
        "Increases quantity by promoting expression"
    )
    UP_REGULATES_QUANTITY_BY_STABILIZATION = (
        "MI:2239",
        "Increases quantity by stabilizing the target"
    )


class ComplexExpansionCv(CvEnum):
    """Complex expansion strategy terms from PSI-MI.

    Describes the method used to expand n-ary complexes into binary interactions
    for analysis and representation purposes.
    """

    # PSI-MI standard complex expansion methods
    BIPARTITE_EXPANSION = (
        "MI:1062",
        "Expand complex into bait-prey binary interactions"
    )
    MATRIX_EXPANSION = (
        "MI:1061",
        "Expand complex into all possible pairwise interactions"
    )
    SPOKE_EXPANSION = (
        "MI:1060",
        "Expand complex into hub-spoke topology"
    )
