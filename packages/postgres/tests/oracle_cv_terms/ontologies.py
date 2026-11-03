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
        "OM:1400",
        "Ontology vocabulary term",
        "Identifies an ontology or controlled vocabulary used by a resource.",
    )

    GENE_ONTOLOGY = ("OM:1401", "Gene Ontology", "https://geneontology.org")
    PSI_MI = ("OM:1402", "Molecular Interactions Ontology", "https://github.com/HUPO-PSI/psi-mi-CV")
    HUMAN_PHENOTYPE_ONTOLOGY = ("OM:1403", "Human Phenotype Ontology", "https://hpo.jax.org/")
    CHEBI = ("OM:1404", "ChEBI", "https://www.ebi.ac.uk/chebi/")
    REACTOME_PATHWAYS = ("OM:1405", "Reactome Pathway Ontology", "https://reactome.org/")
    WIKIPATHWAYS = ("OM:1406", "WikiPathways Ontology", "https://www.wikipathways.org/")
    UNIPROT_KEYWORDS = ("OM:1407", "UniProt Keywords", "https://www.uniprot.org/keywords/")
    OMNIPATH = ("OM:1408", "OmniPath Ontology", "https://omnipathdb.org/")
    MONDO = ("OM:1409", "Mondo Disease Ontology", "https://mondo.monarchinitiative.org/")
    KEGG_PATHWAYS = ("OM:1410", "KEGG Pathway Ontology", "https://www.genome.jp/kegg/pathway.html")
    ENZYME_CLASSIFICATION = ("OM:1411", "Enzyme Classification", "https://www.brenda-enzymes.org/")
    CHEMONT = ("OM:1412", "ChemOnt", "http://classyfire.wishartlab.com/")
    SWISSLIPIDS = ("OM:1413", "SwissLipids", "https://www.swisslipids.org/")
