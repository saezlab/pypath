"""
Parse UniProt data and emit Entity records.

This module converts UniProt protein data into Entity records using the schema
defined in pypath.internals.silver_schema.
"""

from __future__ import annotations

from collections.abc import Generator
import csv

from pypath.share.downloads import download_and_open
from pypath.internals.silver_schema import Entity, Identifier, Annotation
from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    MoleculeAnnotationsCv,
    LicenseCV,
    UpdateCategoryCV,
    ResourceAnnotationCv,
    ResourceCv,
)
from ..internals.tabular_builder import (
    AnnotationsBuilder,
    Column,
    CV,
    EntityBuilder,
    IdentifiersBuilder,
    Map,
)

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


def resource() -> Generator[Entity]:
    """
    Yield resource metadata as an Entity record.

    Yields:
        Entity record with type CV_TERM containing UniProt metadata.
    """
    yield Entity(
        type=EntityTypeCv.CV_TERM,
        identifiers=[
            Identifier(type=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=ResourceCv.UNIPROT),
            Identifier(type=IdentifierNamespaceCv.NAME, value='UniProt'),
        ],
        annotations=[
            Annotation(term=ResourceAnnotationCv.LICENSE, value=str(LicenseCV.CC_BY_4_0)),
            Annotation(term=ResourceAnnotationCv.UPDATE_CATEGORY, value=str(UpdateCategoryCV.REGULAR)),
            Annotation(term=IdentifierNamespaceCv.PUBMED, value='33237286'),
            Annotation(term=ResourceAnnotationCv.URL, value='https://www.uniprot.org/'),
            Annotation(term=ResourceAnnotationCv.DESCRIPTION, value=(
                'UniProt is a comprehensive resource for protein sequence and '
                'functional information. It provides high-quality, manually annotated '
                'protein data including function, structure, localization, interactions, '
                'disease associations, and cross-references to other databases.'
            )),
        ],
    )


def uniprot_proteins() -> Generator[Entity]:
    """
    Download and parse UniProt protein data as Entity records.

    Downloads protein data for human, mouse, and rat from UniProt REST API
    and converts each protein into a Entity with identifiers, annotations,
    membership records, and references.

    Yields:
        Entity records with type PROTEIN
    """
    # Download and open the file
    opener = download_and_open(
        UNIPROT_DATA_URL,
        filename='uniprot_proteins_9606_10090_10116.tsv.gz',
        subfolder='uniprot',
        large=True,
        encoding='utf-8',
        default_mode='r',
        ext='gz',
    )

    # Define the schema mapping
    protein_name_column = Column('Protein names')
    protein_synonym_column = Column('Protein names', delimiter='(')

    schema = EntityBuilder(
        entity_type=EntityTypeCv.PROTEIN,
        identifiers=IdentifiersBuilder(
            # Primary UniProt accession
            CV(term=IdentifierNamespaceCv.UNIPROT, value=Column('Entry')),

            # UniProt entry name
            CV(term=IdentifierNamespaceCv.UNIPROT, value=Column('Entry Name')),

            # Primary gene name
            CV(term=IdentifierNamespaceCv.GENE_NAME_PRIMARY, value=Column('Gene Names (primary)')),

            # Gene synonyms (space-delimited)
            CV(
                term=IdentifierNamespaceCv.GENE_NAME_SYNONYM,
                value=Column('Gene Names (synonym)', delimiter=' '),
            ),

            # Protein name - extract primary (before parentheses)
            CV(
                term=IdentifierNamespaceCv.NAME,
                value=Map(
                    col=protein_name_column,
                    extract=[r'^([^(]+)'],
                ),
            ),

            # Protein name synonyms - split by '(' and extract content before ')'
            # Regex only matches tokens with ')', avoiding the primary name
            CV(
                term=IdentifierNamespaceCv.SYNONYM,
                value=Map(
                    col=protein_synonym_column,
                    extract=[r'^([^)]+)\)'],
                ),
            ),

            # Cross-references (semicolon-delimited)
            CV(term=IdentifierNamespaceCv.ENSEMBL, value=Column('Ensembl', delimiter=';')),
            CV(term=IdentifierNamespaceCv.REFSEQ, value=Column('RefSeq', delimiter=';')),
            CV(term=IdentifierNamespaceCv.PDB, value=Column('PDB', delimiter=';')),
            CV(term=IdentifierNamespaceCv.ALPHAFOLDDB, value=Column('AlphaFoldDB', delimiter=';')),
            CV(term=IdentifierNamespaceCv.KEGG, value=Column('KEGG', delimiter=';')),
            CV(term=IdentifierNamespaceCv.CHEMBL, value=Column('ChEMBL', delimiter=';')),
            CV(term=IdentifierNamespaceCv.SIGNOR, value=Column('SIGNOR', delimiter=';')),
            CV(term=IdentifierNamespaceCv.INTACT, value=Column('IntAct', delimiter=';')),
            CV(term=IdentifierNamespaceCv.BIOGRID, value=Column('BioGRID', delimiter=';')),
            CV(term=IdentifierNamespaceCv.COMPLEXPORTAL, value=Column('ComplexPortal', delimiter=';')),
        ),
        annotations=AnnotationsBuilder(
            # Source annotation
            CV(term=MoleculeAnnotationsCv.SEQUENCE_LENGTH, value=Column('Length')),
            CV(term=MoleculeAnnotationsCv.MASS_DALTON, value=Column('Mass')),

            # Functional annotations
            CV(term=MoleculeAnnotationsCv.FUNCTION, value=Column('Function [CC]')),
            CV(term=MoleculeAnnotationsCv.SUBCELLULAR_LOCATION, value=Column('Subcellular location [CC]')),
            CV(term=MoleculeAnnotationsCv.POST_TRANSLATIONAL_MODIFICATION, value=Column('Post-translational modification')),
            CV(term=MoleculeAnnotationsCv.DISEASE_INVOLVEMENT, value=Column('Involvement in disease')),
            CV(term=MoleculeAnnotationsCv.PATHWAY_PARTICIPATION, value=Column('Pathway')),
            CV(term=MoleculeAnnotationsCv.ACTIVITY_REGULATION, value=Column('Activity regulation')),

            # Feature annotations
            CV(term=MoleculeAnnotationsCv.MUTAGENESIS, value=Column('Mutagenesis')),
            CV(term=MoleculeAnnotationsCv.TRANSMEMBRANE_REGION, value=Column('Transmembrane')),
            CV(term=MoleculeAnnotationsCv.PROTEIN_FAMILY, value=Column('Protein families', delimiter=",")),
            CV(term=MoleculeAnnotationsCv.EC_NUMBER, value=Column('EC number', delimiter=";")),

            # Sequence
            CV(term=MoleculeAnnotationsCv.AMINO_ACID_SEQUENCE, value=Column('Sequence')),

            # GO terms as CV term accessions
            CV(term=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=Column('Gene Ontology IDs', delimiter=';')),

            # Keywords as CV term accessions
            CV(term=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=Column('Keywords IDs', delimiter=';')),

            # Organism as annotation
            CV(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=Column('Organism (ID)')),

            # PubMed references as annotations
            CV(term=IdentifierNamespaceCv.PUBMED, value=Column('PubMed ID', delimiter=';')),
        ),
    )

    # Parse and yield entities
    if opener and opener.result:
        reader = csv.DictReader(opener.result, delimiter='\t')
        for row in reader:
            yield schema(row)
