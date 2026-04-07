"""Ontology controlled vocabulary terms.

This module contains CV terms for ontology vocabularies that may either be
exported directly as ontology artifacts or referenced by annotations in other
resources.
"""

from .core import CvEnum


class OntologyCv(CvEnum):
    """Ontology vocabulary terms used across resources.

    Defines stable OmniPath identifiers for ontology vocabularies referenced by
    resource annotations.
    """

    parent_cv_term = (
        "OM:1200",
        "Ontology vocabulary term",
        "Identifies an ontology or controlled vocabulary used by a resource.",
    )

    GENE_ONTOLOGY = ("OM:1201", "Gene Ontology", "https://geneontology.org")
    PSI_MI = ("OM:1202", "Molecular Interactions Ontology", "https://github.com/HUPO-PSI/psi-mi-CV")
    HUMAN_PHENOTYPE_ONTOLOGY = ("OM:1203", "Human Phenotype Ontology", "https://hpo.jax.org/")
    CHEBI = ("OM:1204", "ChEBI", "https://www.ebi.ac.uk/chebi/")
    REACTOME_PATHWAYS = ("OM:1205", "Reactome Pathway Ontology", "https://reactome.org/")
    WIKIPATHWAYS = ("OM:1206", "WikiPathways Ontology", "https://www.wikipathways.org/")
    UNIPROT_KEYWORDS = ("OM:1207", "UniProt Keywords", "https://www.uniprot.org/keywords/")
    OMNIPATH = ("OM:1208", "OmniPath Ontology", "https://omnipathdb.org/")
