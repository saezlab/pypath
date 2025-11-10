"""Identifier namespace controlled vocabularies.

This module contains all CV terms related to identifier namespaces and reference types,
including database identifiers, chemical structure representations, and reference sources.
"""
from .core import CvEnum


class IdentifierNamespaceCv(CvEnum):
    """Identifier namespace terms backed by PSI-MI and OmniPath accessions.

    Defines standard identifier types used across biological databases.
    Organized by category: proteins/genes, compounds, structures, and names.
    """

    parent_cv_term = "MI:0444"  # database citation - Database citation list names of databases commonly used to cross reference interaction data

    # Protein and gene identifiers (PSI-MI standard)
    UNIPROT = "MI:1097"
    ENTREZ = "MI:0477"
    ENSEMBL = "MI:0476"
    HGNC = "MI:1095"
    REFSEQ = "MI:0481"

    # Compound identifiers (PSI-MI standard)
    CHEBI = "MI:0474"
    PUBCHEM = "MI:0730"
    CHEMBL = "MI:1349"
    CHEMBL_COMPOUND = "MI:0967"
    CHEMBL_TARGET = "MI:1348"
    DRUGBANK = "MI:2002"
    KEGG = "MI:0470"
    KEGG_COMPOUND = "MI:2012"
    CAS = "MI:2011"

    # Structure and interaction databases (OM:0100-0199 range)
    PDB = ("OM:0101", "Protein Data Bank identifier", "https://www.rcsb.org")
    ALPHAFOLDDB = ("OM:0102", "AlphaFold Database identifier", "https://alphafold.ebi.ac.uk")
    INTACT = ("OM:0103", "IntAct interaction database ID", "https://www.ebi.ac.uk/intact")
    BIOGRID = ("OM:0104", "BioGRID interaction database ID", "https://thebiogrid.org")
    COMPLEXPORTAL = ("OM:0105", "Complex Portal identifier", "https://www.ebi.ac.uk/complexportal")
    DIP = ("OM:0106", "Database of Interacting Proteins identifier", "https://dip.doe-mbi.ucla.edu")
    MINT = ("OM:0107", "Molecular INTeraction database identifier", "https://mint.bio.uniroma2.it")
    FLYBASE = ("OM:0108", "FlyBase gene database identifier", "https://flybase.org")
    MIRBASE = ("OM:0109", "miRBase microRNA database identifier", "https://www.mirbase.org")
    RFAM = ("OM:0110", "Rfam RNA families database identifier", "https://rfam.org")
    RNACENTRAL = ("OM:0111", "RNAcentral identifier", "https://rnacentral.org")
    UNIPARC = ("OM:0112", "UniParc (UniProt Archive) identifier", "https://www.uniprot.org/help/uniparc")
    IPI = ("OM:0113", "International Protein Index identifier (deprecated)", "https://www.ebi.ac.uk/proteins/api/doc/")
    IMEX = ("OM:0114", "IMEx consortium identifier", "https://www.imexconsortium.org")
    ENSEMBL_GENOMES = ("OM:0115", "Ensembl Genomes identifier", "https://ensemblgenomes.org")
    GENBANK_NUCL_GI = ("OM:0116", "GenBank nucleotide GI number")
    GENBANK_PROTEIN_GI = ("OM:0117", "GenBank protein GI number")
    GENBANK_IDENTIFIER = ("OM:0118", "GenBank sequence identifier")
    BIND = ("OM:0119", "BIND (Biomolecular Interaction Network Database) identifier", "http://bind.ca")

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
    SMILES = "MI:0239"
    STANDARD_INCHI_KEY = "MI:1101"
    STANDARD_INCHI = "MI:2010"

    # PSI-MI standard reference types
    PUBMED = "MI:0446"
    PUBMED_CENTRAL = "MI:1042"
    DOI = "MI:0574"
    BIORXIV = "MI:2347"