"""Entity type controlled vocabularies.

This module contains CV terms for describing entity types (proteins, genes, compounds, etc.).
Identifier namespaces have been moved to identifiers.py.
"""
from .core import CvEnum


class EntityTypeCv(CvEnum):
    """Top-level entity type classes.

    Only broad PSI-MI and OmniPath entity categories should be listed here.
    More specific classes (protein families, GPCRs, metabolites, etc.) go
    into subordinate subtype vocabularies.
    """

    parent_cv_term = "MI:0313"  # interactor type

    # PSI-MI standard top-level terms
    PROTEIN = "MI:0326"
    GENE = "MI:0250"
    RNA = "MI:0320"
    COMPLEX = "MI:0314"
    SMALL_MOLECULE = "MI:0328"
    PHENOTYPE = "MI:2261"
    STIMULUS = "MI:2260"
    DNA = "MI:0681"
    MACROMOLECULE = "MI:0317"
    NUCLEIC_ACID = "MI:0318"

    # OmniPath high-level types
    PROTEIN_FAMILY = "OM:0010"
    LIPID = ("OM:0011", "Lipid molecule or lipid-like compound")
    CV_TERM = "OM:0012"
    INTERACTION = "OM:0013"
    PATHWAY = "OM:0014"
    REACTION = "OM:0015"
    PHYSICAL_ENTITY = "OM:0016"
    CATALYSIS = "OM:0017"
    CONTROL = "OM:0018"
    DEGRADATION = ("OM:0019", "Degradation reaction")
    FOOD = ("OM:0020", "Food or food product")
    ASSAY = ("OM:0030", "Experimental assay or protocol used to measure biological activity")
    PUBLICATION = ("OM:0031", "A research article, book, patent, or other formal publication.")
    ORGANISM = ("OM:0032", "A living organism, such as a species or strain.")
    CELL_LINE = ("OM:0033", "A cell line or cell culture used in research.")
    TISSUE = ("OM:0034", "A tissue or organ from an organism")


class MoleculeSubtypeCv(CvEnum):
    """Chemical-nature-based molecule subtypes."""

    parent_cv_term = EntityTypeCv.SMALL_MOLECULE

    SYNTHETIC_ORGANIC = ("OM:0029", "Synthetic organic compound")
    NATURAL_PRODUCT = ("OM:0021", "Natural product")
    METABOLITE = ("OM:0022", "Endogenous metabolite")
    INORGANIC = ("OM:0023", "Inorganic compound")
    PEPTIDE = ("OM:0024", "Peptide molecule")
    ANTIBODY = ("OM:0025", "Antibody or immunoglobulin")


class ProteinFunctionalClassCv(CvEnum):
    """Functional subclasses of proteins."""

    parent_cv_term = EntityTypeCv.PROTEIN

    GPCR = ("OM:0040", "G protein-coupled receptor")
    LGIC = ("OM:0041", "Ligand-gated ion channel")
    VGIC = ("OM:0042", "Voltage-gated ion channel")
    OTHER_ION_CHANNEL = ("OM:0043", "Other ion channel")
    ENZYME = ("OM:0044", "Enzyme")
    CATALYTIC_RECEPTOR = ("OM:0045", "Catalytic receptor")
    NUCLEAR_HORMONE_RECEPTOR = ("OM:0046", "Nuclear hormone receptor")
    TRANSPORTER = ("OM:0047", "Transporter protein")
    OTHER_PROTEIN = ("OM:0048", "Other protein")