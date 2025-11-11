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

    # Molecule subtypes by chemical nature (OM:0020-0039 range)
    SYNTHETIC_ORGANIC = ("OM:0020", "Synthetically produced organic compound")
    NATURAL_PRODUCT = ("OM:0021", "Naturally occurring compound isolated from biological sources")
    METABOLITE = ("OM:0022", "Endogenous metabolite or metabolic intermediate")
    INORGANIC = ("OM:0023", "Inorganic compound or ion")
    PEPTIDE = ("OM:0024", "Peptide or small protein molecule")
    ANTIBODY = ("OM:0025", "Antibody or immunoglobulin-based molecule")
    NUCLEIC_ACID = ("OM:0026", "DNA, RNA, or nucleotide-based molecule")

    # Protein functional classes (OM:0040-0059 range)
    GPCR = ("OM:0040", "G protein-coupled receptor")
    LGIC = ("OM:0041", "Ligand-gated ion channel")
    VGIC = ("OM:0042", "Voltage-gated ion channel")
    OTHER_ION_CHANNEL = ("OM:0043", "Other ion channel")
    ENZYME = ("OM:0044", "Enzyme")
    CATALYTIC_RECEPTOR = ("OM:0045", "Catalytic receptor (e.g., receptor tyrosine kinase)")
    NUCLEAR_HORMONE_RECEPTOR = ("OM:0046", "Nuclear hormone receptor")
    TRANSPORTER = ("OM:0047", "Transporter protein")
    OTHER_PROTEIN = ("OM:0048", "Other protein not classified in standard categories")
