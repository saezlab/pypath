"""
Parse UniProt data and emit Entity records.

This module converts UniProt protein data into Entity records using the schema
defined in pypath.internals.silver_schema.
"""

from __future__ import annotations

from collections.abc import Generator

from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    MoleculeAnnotationsCv,
    LicenseCV,
    UpdateCategoryCV,
    ResourceCv,
)
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.base import iter_tsv

# UniProt REST API URL for comprehensive protein data
# Currently hardcoded for human (9606), mouse (10090), and rat (10116)
UNIPROT_DATA_URL = (
    "https://rest.uniprot.org/uniprotkb/stream"
    "?compressed=true"
    "&format=tsv"
    "&query=(taxonomy_id:9606 OR taxonomy_id:10090 OR taxonomy_id:10116) AND reviewed:true"
    "&fields=accession,id,protein_name,length,mass,sequence,gene_primary,gene_synonym,"
    "organism_id,cc_disease,ft_mutagen,cc_subcellular_location,cc_ptm,lit_pubmed_id,"
    "cc_function,xref_ensembl,xref_kegg,cc_pathway,cc_activity_regulation,keywordid,"
    "ec,go_id,ft_transmem,protein_families,xref_refseq,xref_alphafolddb,xref_pdb,"
    "xref_chembl,xref_phosphositeplus,xref_signor,xref_pathwaycommons,xref_intact,"
    "xref_biogrid,xref_complexportal"
)

config = ResourceConfig(
    id=ResourceCv.UNIPROT,
    name='UniProt',
    url='https://www.uniprot.org/',
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='33237286',
    description=(
        'UniProt is a comprehensive resource for protein sequence and '
        'functional information. It provides high-quality, manually annotated '
        'protein data including function, structure, localization, interactions, '
        'disease associations, and cross-references to other databases.'
    ),
)

f = FieldConfig(
    extract={
        'protein_name': r'^([^(]+)',
        'protein_synonym': r'^([^)]+)\)',
    },
)

protein_name_column = f('Protein names', extract='protein_name')
protein_synonym_column = f('Protein names', delimiter='(', extract='protein_synonym')

proteins_schema = EntityBuilder(
    entity_type=EntityTypeCv.PROTEIN,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.UNIPROT, value=f('Entry')),
        CV(term=IdentifierNamespaceCv.UNIPROT, value=f('Entry Name')),
        CV(term=IdentifierNamespaceCv.GENE_NAME_PRIMARY, value=f('Gene Names (primary)')),
        CV(
            term=IdentifierNamespaceCv.GENE_NAME_SYNONYM,
            value=f('Gene Names (synonym)', delimiter=' '),
        ),
        CV(
            term=IdentifierNamespaceCv.NAME,
            value=protein_name_column,
        ),
        CV(
            term=IdentifierNamespaceCv.SYNONYM,
            value=protein_synonym_column,
        ),
        CV(term=IdentifierNamespaceCv.ENSEMBL, value=f('Ensembl', delimiter=';')),
        CV(term=IdentifierNamespaceCv.REFSEQ, value=f('RefSeq', delimiter=';')),
        CV(term=IdentifierNamespaceCv.PDB, value=f('PDB', delimiter=';')),
        CV(term=IdentifierNamespaceCv.ALPHAFOLDDB, value=f('AlphaFoldDB', delimiter=';')),
        CV(term=IdentifierNamespaceCv.KEGG, value=f('KEGG', delimiter=';')),
        CV(term=IdentifierNamespaceCv.CHEMBL, value=f('ChEMBL', delimiter=';')),
        CV(term=IdentifierNamespaceCv.SIGNOR, value=f('SIGNOR', delimiter=';')),
        CV(term=IdentifierNamespaceCv.INTACT, value=f('IntAct', delimiter=';')),
        CV(term=IdentifierNamespaceCv.BIOGRID, value=f('BioGRID', delimiter=';')),
        CV(term=IdentifierNamespaceCv.COMPLEXPORTAL, value=f('ComplexPortal', delimiter=';')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=MoleculeAnnotationsCv.SEQUENCE_LENGTH, value=f('Length')),
        CV(term=MoleculeAnnotationsCv.MASS_DALTON, value=f('Mass')),
        CV(term=MoleculeAnnotationsCv.FUNCTION, value=f('Function [CC]')),
        CV(term=MoleculeAnnotationsCv.SUBCELLULAR_LOCATION, value=f('Subcellular location [CC]')),
        CV(term=MoleculeAnnotationsCv.POST_TRANSLATIONAL_MODIFICATION, value=f('Post-translational modification')),
        CV(term=MoleculeAnnotationsCv.DISEASE_INVOLVEMENT, value=f('Involvement in disease')),
        CV(term=MoleculeAnnotationsCv.PATHWAY_PARTICIPATION, value=f('Pathway')),
        CV(term=MoleculeAnnotationsCv.ACTIVITY_REGULATION, value=f('Activity regulation')),
        CV(term=MoleculeAnnotationsCv.MUTAGENESIS, value=f('Mutagenesis')),
        CV(term=MoleculeAnnotationsCv.TRANSMEMBRANE_REGION, value=f('Transmembrane')),
        CV(term=MoleculeAnnotationsCv.PROTEIN_FAMILY, value=f('Protein families', delimiter=',')),
        CV(term=MoleculeAnnotationsCv.EC_NUMBER, value=f('EC number', delimiter=';')),
        CV(term=MoleculeAnnotationsCv.AMINO_ACID_SEQUENCE, value=f('Sequence')),
        CV(term=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=f('Gene Ontology IDs', delimiter=';')),
        CV(term=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=f('Keywords IDs', delimiter=';')),
        CV(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=f('Organism (ID)')),
        CV(term=IdentifierNamespaceCv.PUBMED, value=f('PubMed ID', delimiter=';')),
    ),
)

resource = Resource(
    config,
    proteins=Dataset(
        download=Download(
            url=UNIPROT_DATA_URL,
            filename='uniprot_proteins_9606_10090_10116.tsv.gz',
            subfolder='uniprot',
            large=True,
            encoding='utf-8',
            default_mode='r',
            ext='gz',
        ),
        mapper=proteins_schema,
        raw_parser=iter_tsv,
    ),
)
