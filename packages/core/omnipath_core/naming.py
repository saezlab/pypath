"""Naming and namespace normalization utilities for OmniPath."""

from __future__ import annotations

from enum import Enum, StrEnum
from typing import Any


class Namespace(StrEnum):
    """Canonical identifier namespaces for OmniPath.

    Subclasses StrEnum so all members are directly instances of str,
    compatible with PyArrow schemas, dictionary keys, and string formatting,
    while offering IDE autocompletion and static type checking.
    """

    # Proteins and Genes
    UNIPROT = "uniprot"
    UNIPROT_TREMBL = "uniprot_trembl"
    UNIPROT_SEC = "uniprot-sec"
    UNIPROT_ENTRY = "uniprot_entry"
    GENESYMBOL = "genesymbol"
    GENESYMBOL_SYN = "genesymbol-syn"
    ENTREZ = "entrez"
    ENSG = "ensg"
    ENST = "enst"
    ENSP = "ensp"
    HGNC = "hgnc"
    REFSEQ = "refseq"
    REFSEQ_PROTEIN = "refseq_protein"

    # Small Molecules & Chemicals
    CHEBI = "chebi"
    PUBCHEM = "pubchem"
    PUBCHEM_SUBSTANCE = "pubchem_substance"
    CHEMBL = "chembl"
    CHEMBL_TARGET = "chembl_target"
    CHEMBL_ACTIVITY = "chembl_activity"
    CHEMBL_MECHANISM = "chembl_mechanism"
    DRUGBANK = "drugbank"
    KEGG = "kegg"
    KEGG_REACTION = "kegg_reaction"
    KEGG_PATHWAY = "kegg_pathway"
    CAS = "cas"
    HMDB = "hmdb"
    LIPIDMAPS = "lipidmaps"
    SWISSLIPIDS = "swisslipids"
    METANETX = "metanetx"
    BINDINGDB = "bindingdb"
    GUIDETOPHARMA = "guidetopharma"
    ZINC = "zinc"
    DRUGCENTRAL = "drugcentral"
    INCHIKEY = "inchikey"
    INCHI = "inchi"
    SMILES = "smiles"

    # Complexes & Interactions
    SIGNOR = "signor"
    COMPLEXPORTAL = "complexportal"
    INTACT = "intact"
    CORUM = "corum"
    BIOGRID = "biogrid"
    DIP = "dip"
    MINT = "mint"
    STRING = "string"

    # Structures
    PDB = "pdb"
    ALPHAFOLDDB = "alphafolddb"

    # RNA & Non-coding
    MIRBASE = "mirbase"
    MIRBASE_PRECURSOR = "mirbase_precursor"
    MIRBASE_MATURE = "mirbase_mature"
    RFAM = "rfam"
    RNACENTRAL = "rnacentral"
    UNIPARC = "uniparc"

    # Ontologies & Pathways
    GO = "go"
    UNIPROT_KEYWORD = "uniprot_keyword"
    REACTOME = "reactome"
    WIKIPATHWAYS = "wikipathways"
    DOID = "doid"
    RAMP = "ramp"
    EC = "ec"
    PFOCR = "pfocr"
    REFMET = "refmet"
    RHEA = "rhea"
    ECOCYC = "ecocyc"
    METACYC = "metacyc"
    TCDB = "tcdb"

    # Literature & Annotations
    PUBMED = "pubmed"
    PUBMED_CENTRAL = "pubmed_central"
    DOI = "doi"
    BIORXIV = "biorxiv"

    # Source-owned identifier systems (distinct accession domains stay separate).
    FOODB = "foodb"
    PHENOL_EXPLORER_FOOD = "phenol_explorer_food"
    PHENOL_EXPLORER_COMPOUND = "phenol_explorer_compound"
    PTFI = "ptfi"
    FOODON = "foodon"
    HUMAN_GEM_METABOLITE = "human_gem_metabolite"
    HUMAN_GEM_REACTION = "human_gem_reaction"
    CHEMONT = "chemont"
    CHEMSPIDER = "chemspider"
    CELLCHAT = "cellchat"
    CELLINKER = "cellinker"
    CDB = "cdb"
    MEBOCOST = "mebocost"
    MRCLINKSDB = "mrclinksdb"
    BIGG_METABOLITE = "bigg_metabolite"
    ENVIPATH = "envipath"
    KEGG_DRUG = "kegg_drug"
    KEGG_GLYCAN = "kegg_glycan"
    SABIORK_COMPOUND = "sabiork_compound"
    SEED_COMPOUND = "seed_compound"
    GUIDETOPHARMA_TARGET = "guidetopharma_target"
    BIGG_REACTION = "bigg_reaction"
    METANETX_REACTION = "metanetx_reaction"
    REACTOME_ID = "reactome_id"
    KEGG_ORTHOLOGY = "kegg_orthology"
    KEGG_GENE = "kegg_gene"
    KEGG_RCLASS = "kegg_rclass"
    WIKIPATHWAYS_VERSION = "wikipathways_version"
    BIND = "bind"
    GENBANK = "genbank"
    GENBANK_GI = "genbank_gi"
    ENSEMBL = "ensembl"
    ENSEMBL_GENOMES = "ensembl_genomes"
    FLYBASE = "flybase"
    IMEX = "imex"
    IPI = "ipi"
    MI = "mi"
    CHEMBL_INTERNAL_ID = "chembl_internal_id"
    MACDB_TRAIT = "macdb_trait"
    HPO = "hpo"
    MONDO = "mondo"
    BTO = "bto"
    CL = "cl"
    UBERON = "uberon"
    APO = "apo"
    PATO = "pato"
    MP = "mp"
    WBPHENOTYPE = "wbphenotype"
    ECO = "eco"
    INTERPRO = "interpro"
    PFAM = "pfam"
    UMLS = "umls"
    ONCOTREE = "oncotree"
    MGI = "mgi"
    NCIT = "ncit"
    RXCUI = "rxcui"
    WIKIDATA = "wikidata"
    HEMONC = "hemonc"
    DRUGSATFDA_NDA = "drugsatfda_nda"
    DRUGSATFDA_ANDA = "drugsatfda_anda"
    CHEMIDPLUS = "chemidplus"
    MSIGDB = "msigdb"
    OM = "om"

    # Generic & Entity Metadata
    NAME = "name"
    SYNONYM = "synonym"
    NCBI_TAX_ID = "ncbi_tax_id"
    CV_TERM = "cv_term"


def normalize_namespace(namespace: Any) -> Namespace | str:
    """Read a namespace without enum aliases or identifier-content inference.

    Unregistered source prefixes are preserved for evidence. Add a Namespace member
    when a source becomes supported; names/taxonomy/publications use Biolink slots.
    """
    if type(namespace) is str:
        raw = namespace.strip().lower()
        return Namespace._value2member_map_.get(raw, raw)
    if namespace is None:
        return ""
    if isinstance(namespace, Namespace):
        return namespace
    if isinstance(namespace, Enum) or not isinstance(namespace, str):
        raise ValueError(f"Expected Namespace or its serialized string, got {namespace!r}")
    raw = namespace.strip().lower()
    return Namespace._value2member_map_.get(raw, raw)
