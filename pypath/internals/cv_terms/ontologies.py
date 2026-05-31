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
        "OM:1240",
        "Ontology vocabulary term",
        "Identifies an ontology or controlled vocabulary used by a resource.",
    )

    GENE_ONTOLOGY = ("OM:1241", "Gene Ontology", "https://geneontology.org")
    PSI_MI = ("OM:1242", "Molecular Interactions Ontology", "https://github.com/HUPO-PSI/psi-mi-CV")
    HUMAN_PHENOTYPE_ONTOLOGY = ("OM:1243", "Human Phenotype Ontology", "https://hpo.jax.org/")
    CHEBI = ("OM:1244", "ChEBI", "https://www.ebi.ac.uk/chebi/")
    REACTOME_PATHWAYS = ("OM:1245", "Reactome Pathway Ontology", "https://reactome.org/")
    WIKIPATHWAYS = ("OM:1246", "WikiPathways Ontology", "https://www.wikipathways.org/")
    UNIPROT_KEYWORDS = ("OM:1247", "UniProt Keywords", "https://www.uniprot.org/keywords/")
    OMNIPATH = ("OM:1248", "OmniPath Ontology", "https://omnipathdb.org/")
    MONDO = ("OM:1249", "Mondo Disease Ontology", "https://mondo.monarchinitiative.org/")
    KEGG_PATHWAYS = ("OM:1250", "KEGG Pathway Ontology", "https://www.genome.jp/kegg/pathway.html")
    ENZYME_CLASSIFICATION = ("OM:1251", "Enzyme Classification", "https://www.brenda-enzymes.org/")
    CHEMONT = ("OM:1252", "ChemOnt", "http://classyfire.wishartlab.com/")
