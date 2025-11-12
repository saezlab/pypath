"""Resource controlled vocabularies.

This module contains CV terms for biological data sources and resources.
Each resource is represented as a controlled vocabulary term that can be used
to identify and describe data sources in the OmniPath ecosystem.
"""
from .core import CvEnum


class ResourceCv(CvEnum):
    """Resource/data source controlled vocabulary terms.

    Defines standard identifiers for biological databases and data sources.
    Each resource has a unique OmniPath accession in the OM:1100-1999 range.
    """

    parent_cv_term = ("OM:1100", "Resource term", "Identifies a biological database or data source.")

    # Ontology resources (OM:1100-1119 range)
    OMNIPATH_ONTOLOGY = ("OM:1100", "OmniPath Ontology", "https://omnipathdb.org/")
    UNIPROT_KEYWORDS = ("OM:1101", "UniProt Keywords", "https://www.uniprot.org/keywords/")
    GENE_ONTOLOGY = ("OM:1102", "Gene Ontology", "https://geneontology.org")
    PSI_MI = ("OM:1103", "PSI-MI Controlled Vocabulary", "https://github.com/HUPO-PSI/psi-mi-CV")

    # Protein databases (OM:1120-1149 range)
    UNIPROT = ("OM:1120", "UniProt", "https://www.uniprot.org")

    # Interaction databases (OM:1150-1179 range)
    INTACT = ("OM:1150", "IntAct", "https://www.ebi.ac.uk/intact")
    SIGNOR = ("OM:1152", "SIGNOR", "https://signor.uniroma2.it/")
    BINDINGDB = ("OM:1153", "BindingDB", "https://www.bindingdb.org/")

    # Metabolite and lipid databases (OM:1170-1179 range)
    HMDB = ("OM:1170", "Human Metabolome Database", "https://hmdb.ca/")
    LIPIDMAPS = ("OM:1171", "LIPID MAPS Structure Database", "https://lipidmaps.org/")
    SWISSLIPIDS = ("OM:1172", "SwissLipids", "https://www.swisslipids.org/")

    # Pharmacological databases (OM:1180-1199 range)
    GUIDETOPHARMA = ("OM:1180", "Guide to Pharmacology", "https://www.guidetopharmacology.org")
