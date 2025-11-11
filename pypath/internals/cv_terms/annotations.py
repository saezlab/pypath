"""Annotation-related controlled vocabularies.

This module contains all CV terms used for annotating interactions and entities,
including roles, interaction types, detection methods, biological effects, causal
mechanisms and statements, complex expansion methods, and curation metadata.
"""
from .core import CvEnum


class MembershipRoleCv(CvEnum):
    """Relationship types for entity membership (is_member_of field).

    Defines how entities relate to their parent entities in hierarchical
    or compositional relationships (e.g., proteins in complexes, ontology terms).
    """

    parent_cv_term = ("OM:0300", "Membership role term", "Describes the membership relationship of an entity.")  # OmniPath-specific term - no PSI-MI parent

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


class PharmacologicalActionCv(CvEnum):
    """Pharmacological action terms describing drug-target interactions.

    Describes the pharmacological effect or mechanism of action of a ligand
    on its target, commonly used in drug databases like Guide to Pharmacology.
    These terms describe both the functional outcome and mechanism.
    """

    parent_cv_term = ("OM:0900", "Pharmacological action term", "Describes the pharmacological effect or mechanism of a ligand on its target.")

    # Agonist actions (OM:0901-0919 range)
    AGONIST = ("OM:0901", "Ligand that activates a receptor")
    FULL_AGONIST = ("OM:0902", "Ligand that produces maximal receptor activation")
    PARTIAL_AGONIST = ("OM:0903", "Ligand that produces submaximal receptor activation")
    INVERSE_AGONIST = ("OM:0904", "Ligand that reduces constitutive receptor activity")
    BIASED_AGONIST = ("OM:0905", "Ligand that selectively activates specific signaling pathways")
    IRREVERSIBLE_AGONIST = ("OM:0906", "Agonist that forms permanent or very stable binding")

    # Antagonist actions (OM:0920-0929 range)
    ANTAGONIST = ("OM:0920", "Ligand that blocks receptor activation")
    COMPETITIVE = ("OM:0921", "Antagonist that competes with agonist for binding site")
    NON_COMPETITIVE = ("OM:0922", "Antagonist that binds to different site than agonist")

    # Activation/Inhibition (OM:0930-0949 range)
    ACTIVATION = ("OM:0930", "General activation of target activity")
    INHIBITION = ("OM:0931", "General inhibition of target activity")
    IRREVERSIBLE_INHIBITION = ("OM:0932", "Permanent or very stable inhibition")
    FEEDBACK_INHIBITION = ("OM:0933", "Inhibition through feedback mechanism")

    # Modulator actions (OM:0950-0969 range)
    POSITIVE = ("OM:0950", "Positive modulation of target activity")
    NEGATIVE = ("OM:0951", "Negative modulation of target activity")
    POTENTIATION = ("OM:0952", "Enhancement of target response")
    NEUTRAL = ("OM:0953", "Binding without functional effect")

    # Channel-specific actions (OM:0970-0979 range)
    PORE_BLOCKER = ("OM:0970", "Blocks ion channel pore")
    SLOWS_INACTIVATION = ("OM:0971", "Delays channel inactivation")
    VOLTAGE_DEPENDENT_INHIBITION = ("OM:0972", "Inhibition dependent on membrane voltage")

    # Other/mixed (OM:0980-0999 range)
    BINDING = ("OM:0980", "Simple binding without specified functional outcome")
    BIPHASIC = ("OM:0981", "Dual or concentration-dependent effects")
    MIXED = ("OM:0982", "Multiple or mixed pharmacological actions")
    UNKNOWN = ("OM:0983", "Pharmacological action not determined")
    NONE = ("OM:0984", "No pharmacological action")


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


class CurationCv(CvEnum):
    """Curation and annotation metadata terms.

    Terms for capturing metadata about the curation process and
    supporting evidence for interactions and entities.
    """

    parent_cv_term = ("OM:0400", "Curation term", "Describes the results of a curation process and supporting evidence for interactions and entities.")  # OmniPath-specific term - no PSI-MI parent

    # OmniPath curation terms (OM:0400-0499 range)
    EVIDENCE_SENTENCE = (
        "OM:0401",
        "Sentence or text excerpt from literature supporting the interaction or annotation"
    )

    # PSI-MI standard curation terms
    COMMENT = "MI:0612"

class MoleculeAnnotationsCv(CvEnum):
    """Controlled vocabulary for molecule annotation types.

    Describes various annotation categories used in molecule entries,
    such as domains, families, transmembrane regions, and mutagenesis data.
    """
    parent_cv_term = ("OM:0600", "Molecule annotation type term", "Describes various annotation categories used in molecule entries.")
    AMINO_ACID_SEQUENCE = ("OM:0600", "Sequence of amino acids of the protein or peptide")
    SEQUENCE_LENGTH = ("OM:0601", "Length of the protein sequence")
    MASS_DALTON = ("OM:0602", "Molecular mass in Daltons")
    FUNCTION = ("OM:0603", "Functional description of the molecule")
    SUBCELLULAR_LOCATION = ("OM:0604", "Subcellular localization information")
    POST_TRANSLATIONAL_MODIFICATION = ("OM:0605", "Post-translational modification details")
    DISEASE_INVOLVEMENT = ("OM:0606", "Information about involvement in diseases")
    PATHWAY_PARTICIPATION = ("OM:0607", "Pathways in which the molecule is involved")
    ACTIVITY_REGULATION = ("OM:0608", "Regulation of molecular activity")
    TRANSMEMBRANE_REGION = ("OM:0609", "Information about transmembrane regions")
    PROTEIN_FAMILY = ("OM:0610", "Protein family classification")
    EC_NUMBER = ("OM:0611", "Enzyme Commission number")
    MUTAGENESIS = ("OM:0612", "Details about mutagenesis experiments")
    DESCRIPTION = ("OM:0613", "General description of the molecule")

    # Lipid classification terms (OM:0614-0622 range)
    LIPID_CATEGORY = ("OM:0614", "Lipid category classification (e.g., Fatty Acyls, Glycerolipids)")
    LIPID_MAIN_CLASS = ("OM:0615", "Main lipid class within a category")
    LIPID_SUB_CLASS = ("OM:0616", "Sub-class within a lipid main class")
    LIPID_HIERARCHY_LEVEL = ("OM:0619", "Level of specificity in lipid classification hierarchy (e.g., Species, Isomer, Class)")
    LIPID_STRUCTURAL_COMPONENTS = ("OM:0620", "Structural components of the lipid molecule (e.g., fatty acid chains)")

    # Molecular properties (OM:0623-0629 range)
    MOLECULAR_CHARGE = ("OM:0623", "Electric charge of the molecule at specified pH")
    PRIMARY_TARGET = ("OM:0624", "Indicates whether this is the primary/main target of a ligand (true/false)")
    ENDOGENOUS = ("OM:0625", "Indicates whether the ligand is an endogenous molecule for the target organism (true/false)")

    # Affinity measurements (OM:0626-0628 range) - used with AffinityUnitCv
    AFFINITY_HIGH = ("OM:0626", "High affinity measurement value")
    AFFINITY_LOW = ("OM:0627", "Low affinity measurement value")
    AFFINITY_MEDIAN = ("OM:0628", "Median affinity measurement value")

    # Drug/molecule status and properties (OM:0650-0654 range)
    APPROVED = ("OM:0650", "Indicates whether the drug is approved for clinical use")
    WITHDRAWN = ("OM:0651", "Indicates whether the drug has been withdrawn from market")
    LABELLED = ("OM:0652", "Indicates whether the molecule is isotopically or fluorescently labelled")
    RADIOACTIVE = ("OM:0653", "Indicates whether the molecule contains radioactive isotopes")
    ANTIBACTERIAL = ("OM:0654", "Indicates whether the molecule has antibacterial activity")


class InteractionParameterCv(CvEnum):
    """Interaction parameter terms from PSI-MI.

    Describes kinetic and thermodynamic parameters for enzymatic or binding
    studies, including affinity measurements, rate constants, and experimental
    conditions (e.g., pH, temperature).
    """

    parent_cv_term = "MI:0640"  # parameter type - Parameter for enzymatic or binding kinetic studies

    # Binding affinity measurements (equilibrium constants)
    KI = "MI:0643"  # Equilibrium constant for dissociation of an inhibitor. Unit Molar.
    KD = "MI:0646"  # The equilibrium dissociation constant. Unit Molar.
    IC50 = "MI:0641"  # Molar concentration producing 50% of maximum inhibitory response. Unit Molar.
    EC50 = "MI:0642"  # Molar concentration producing 50% of maximum response for agonist. Unit Molar.

    # Rate constants
    KON = "MI:0834"  # Association rate constant or rate of complex formation. Unit M-1 s-1
    KOFF = "MI:0835"  # Dissociation rate constant measuring stability of a complex. Unit s-1

    # Experimental conditions
    PH = "MI:0837"  # pH at which interaction was determined
    TEMPERATURE = "MI:0836"  # Temperature at which interaction was determined. Unit KELVIN (K)
    TEMPERATURE_CELSIUS = ("OM:0701", "Temperature at which interaction was determined. Unit CELSIUS (C)")


class AffinityUnitCv(CvEnum):
    """Affinity measurement unit terms.

    Describes the units used for expressing binding affinity measurements.
    Most commonly logarithmic scales (pKd, pIC50, etc.) or direct measurements.
    """

    parent_cv_term = ("OM:0700", "Affinity unit term", "Describes the unit of measurement for binding affinity values.")

    # Logarithmic affinity units (OM:0702-0719 range)
    PKD = ("OM:0702", "Negative logarithm of dissociation constant (pKd = -log10(Kd))")
    PKI = ("OM:0703", "Negative logarithm of inhibition constant (pKi = -log10(Ki))")
    PIC50 = ("OM:0704", "Negative logarithm of IC50 (pIC50 = -log10(IC50))")
    PEC50 = ("OM:0705", "Negative logarithm of EC50 (pEC50 = -log10(EC50))")
    PKB = ("OM:0706", "Negative logarithm of equilibrium constant for antagonist (pKb = -log10(Kb))")
    PA2 = ("OM:0707", "Negative logarithm of molar concentration of antagonist that makes it necessary to double agonist concentration")

    # Direct measurement units would go here if needed (OM:0720-0739 range)
    # MOLAR = ("OM:0720", "Molar concentration (M)")
    # NANOMOLAR = ("OM:0721", "Nanomolar concentration (nM)")


class LigandTypeCv(CvEnum):
    """Ligand type classification terms.

    Describes the molecular type or classification of a ligand based on
    its structure, origin, or functional characteristics.
    """

    parent_cv_term = ("OM:1000", "Ligand type term", "Describes the molecular type or classification of a ligand.")

    # Functional types (OM:1001-1019 range)
    AGONIST = ("OM:1001", "Ligand that activates a receptor")
    ANTAGONIST = ("OM:1002", "Ligand that blocks receptor activation")
    ACTIVATOR = ("OM:1003", "Molecule that increases target activity")
    INHIBITOR = ("OM:1004", "Molecule that decreases target activity")
    ALLOSTERIC_MODULATOR = ("OM:1005", "Ligand that binds at site distinct from active site and modulates activity")

    # Channel-specific types (OM:1020-1029 range)
    CHANNEL_BLOCKER = ("OM:1020", "Molecule that blocks ion channel activity")
    GATING_INHIBITOR = ("OM:1021", "Molecule that inhibits channel gating mechanism")

    # Structural/molecular types (OM:1030-1049 range)
    ANTIBODY = ("OM:1030", "Immunoglobulin-based ligand")
    FUSION_PROTEIN = ("OM:1031", "Chimeric protein combining domains from different sources")

    # Specificity types (OM:1050-1059 range)
    SUBUNIT_SPECIFIC = ("OM:1050", "Ligand that selectively targets specific subunit(s)")

    # Unspecified (OM:1099)
    NONE = ("OM:1099", "No specific ligand type classification")


class OntologyAnnotationCv(CvEnum):
    """Ontology annotation terms.

    Describes metadata for ontology terms (CV_TERM entities).
    Note: Relationships like is_a and xref should use CV_TERM_ACCESSION in identifiers.
    """

    parent_cv_term = ("OM:0800", "Ontology annotation term", "Describes metadata for ontology terms.")

    # Core ontology metadata (OM:0800-0899 range)
    DEFINITION = ("OM:0801", "Textual definition of an ontology term")
    COMMENT = ("OM:0805", "Additional comment or note about an ontology term")
    IS_OBSOLETE = ("OM:0806", "Indicates whether the term is obsolete")
