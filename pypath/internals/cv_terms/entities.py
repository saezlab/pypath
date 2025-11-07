"""Entity type and identifier namespace controlled vocabularies.

This module contains CV terms for describing entities (proteins, genes, compounds, etc.)
and their identifier namespaces (UniProt, ChEBI, PubChem, etc.).
"""
from .core import CvEnum


class EntityTypeCv(CvEnum):
    """Common PSI-MI and OmniPath entity type terms.

    Defines the types of entities that can be represented in the system,
    including standard PSI-MI types and OmniPath-specific extensions.
    """

    # PSI-MI standard terms
    PROTEIN = ("MI:0326", "Protein molecule or polypeptide chain")
    GENE = ("MI:0250", "Gene entity or genetic locus")
    RNA = ("MI:0320", "RNA molecule including mRNA, ncRNA, etc.")
    PROTEIN_COMPLEX = ("MI:0315", "Complex of multiple protein subunits")
    SMALL_MOLECULE = ("MI:0328", "Small chemical compound or drug")
    PHENOTYPE = ("MI:2261", "Phenotypic effect or observable trait")
    STIMULUS = ("MI:2260", "External stimulus or experimental condition")

    # OmniPath-specific terms (OM:0010-0099 range)
    PROTEIN_FAMILY = ("OM:0010", "Protein family grouping related proteins")
    LIPID = ("OM:0011", "Lipid molecule or lipid-like compound")
    CV_TERM = ("OM:0012", "Controlled vocabulary term entity")
    INTERACTION = ("OM:0013", "Interaction entity for representing interactions as nodes")


class IdentifierNamespaceCv(CvEnum):
    """Identifier namespace terms backed by PSI-MI and OmniPath accessions.

    Defines standard identifier types used across biological databases.
    Organized by category: proteins/genes, compounds, structures, and names.
    """

    # Protein and gene identifiers (PSI-MI standard)
    UNIPROT = ("MI:1097", "UniProt protein accession", "https://www.uniprot.org")
    ENTREZ = ("MI:0477", "NCBI Entrez Gene ID", "https://www.ncbi.nlm.nih.gov/gene")
    ENSEMBL = ("MI:0476", "Ensembl gene/transcript/protein ID", "https://www.ensembl.org")
    HGNC = ("MI:1095", "HUGO Gene Nomenclature Committee ID", "https://www.genenames.org")
    REFSEQ = ("MI:0481", "NCBI RefSeq accession", "https://www.ncbi.nlm.nih.gov/refseq")

    # Compound identifiers (PSI-MI standard)
    CHEBI = ("MI:0474", "Chemical Entities of Biological Interest ID", "https://www.ebi.ac.uk/chebi")
    PUBCHEM = ("MI:0730", "PubChem database identifier", "https://pubchem.ncbi.nlm.nih.gov")
    CHEMBL = ("MI:1349", "ChEMBL database identifier", "https://www.ebi.ac.uk/chembl")
    CHEMBL_COMPOUND = ("MI:0967", "ChEMBL compound ID", "https://www.ebi.ac.uk/chembl")
    CHEMBL_TARGET = ("MI:1348", "ChEMBL target ID", "https://www.ebi.ac.uk/chembl")
    DRUGBANK = ("MI:2002", "DrugBank database identifier", "https://www.drugbank.ca")
    KEGG = ("MI:0470", "KEGG database identifier", "https://www.genome.jp/kegg")
    KEGG_COMPOUND = ("MI:2012", "KEGG Compound ID", "https://www.genome.jp/kegg/compound")
    CAS = ("MI:2011", "Chemical Abstracts Service Registry Number", "https://www.cas.org")

    # Structure and interaction databases (OM:0100-0199 range)
    PDB = ("OM:0101", "Protein Data Bank identifier", "https://www.rcsb.org")
    ALPHAFOLDDB = ("OM:0102", "AlphaFold Database identifier", "https://alphafold.ebi.ac.uk")
    INTACT = ("OM:0103", "IntAct interaction database ID", "https://www.ebi.ac.uk/intact")
    BIOGRID = ("OM:0104", "BioGRID interaction database ID", "https://thebiogrid.org")
    COMPLEXPORTAL = ("OM:0105", "Complex Portal identifier", "https://www.ebi.ac.uk/complexportal")

    # Additional specialized identifiers (OM:0001-0009 range)
    REFSEQ_PROTEIN = ("OM:0001", "NCBI RefSeq protein accession (NP_* format)")
    PUBCHEM_COMPOUND = ("OM:0002", "PubChem Compound ID (CID)")
    LIPIDMAPS = ("OM:0003", "LIPID MAPS structure database ID", "https://www.lipidmaps.org")
    HMDB = ("OM:0004", "Human Metabolome Database ID", "https://hmdb.ca")
    METANETX = ("OM:0005", "MetaNetX/MNXref chemical ID", "https://www.metanetx.org")
    BINDINGDB = ("OM:0006", "BindingDB database identifier", "https://www.bindingdb.org")
    SIGNOR = ("OM:0007", "SIGNOR database identifier", "https://signor.uniroma2.it")
    GUIDETOPHARMA = ("OM:0008", "Guide to Pharmacology identifier", "https://www.guidetopharmacology.org")
    SWISSLIPIDS = ("OM:0009", "SwissLipids database ID", "https://www.swisslipids.org")

    # Gene and protein names (OM:0200-0209 range)
    GENE_NAME_PRIMARY = ("OM:0200", "Primary gene name or symbol")
    GENE_NAME_SYNONYM = ("OM:0201", "Alternative gene name or synonym")
    NAME = ("OM:0202", "Generic primary name for any entity type")
    SYNONYM = ("OM:0203", "Generic synonym or alternative name")
    CV_TERM_ACCESSION = ("OM:0204", "Controlled vocabulary term accession")
    NCBI_TAX_ID = ("OM:0205", "NCBI Taxonomy database ID", "https://www.ncbi.nlm.nih.gov/taxonomy")

    # Structure representations (PSI-MI standard)
    SMILES = ("MI:0239", "Simplified Molecular Input Line Entry System")
    STANDARD_INCHI_KEY = ("MI:1101", "IUPAC International Chemical Identifier InChIKey")
    STANDARD_INCHI = ("MI:2010", "IUPAC International Chemical Identifier InChI")
