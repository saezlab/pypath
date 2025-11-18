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
    PROTEIN_COMPLEX = "MI:0315"
    SMALL_MOLECULE = "MI:0328"
    PHENOTYPE = "MI:2261"
    STIMULUS = "MI:2260"

    # OmniPath high-level types
    PROTEIN_FAMILY = "OM:0010"
    LIPID = ("OM:0011", "Lipid molecule or lipid-like compound")
    CV_TERM = "OM:0012"
    INTERACTION = "OM:0013"
    PATHWAY = "OM:0014"
    REACTION = "OM:0015"

class MoleculeSubtypeCv(CvEnum):
    """Chemical-nature-based molecule subtypes."""

    parent_cv_term = EntityTypeCv.SMALL_MOLECULE

    SYNTHETIC_ORGANIC = ("OM:0020", "Synthetic organic compound")
    NATURAL_PRODUCT = ("OM:0021", "Natural product")
    METABOLITE = ("OM:0022", "Endogenous metabolite")
    INORGANIC = ("OM:0023", "Inorganic compound")
    PEPTIDE = ("OM:0024", "Peptide molecule")
    ANTIBODY = ("OM:0025", "Antibody or immunoglobulin")
    NUCLEIC_ACID = ("OM:0026", "DNA/RNA/nucleotide-based molecule")


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