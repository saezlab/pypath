#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright 2014-2025
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: see the file `README.rst`
#  Contact: Dénes Türei (turei.denes@gmail.com)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      https://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: https://pypath.omnipathdb.org/
#

"""
Parse UniProt data and emit Entity records.

This module converts UniProt protein data into Entity records using the schema
defined in pypath.internals.silver_schema.
"""

from __future__ import annotations

from collections.abc import Generator
import csv

from pypath.share.downloads import download_and_open
from pypath.internals.silver_schema import Entity as SilverEntity, Identifier, Annotation
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
    Annotations,
    Column,
    Entity,
    Identifiers,
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


def uniprot() -> Generator[SilverEntity]:
    """
    Yield resource metadata as an Entity record.

    Yields:
        Entity record with type CV_TERM containing UniProt metadata.
    """
    yield SilverEntity(
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


def uniprot_proteins() -> Generator[SilverEntity]:
    """
    Download and parse UniProt protein data as Entity records.

    Downloads protein data for human, mouse, and rat from UniProt REST API
    and converts each protein into a SilverEntity with identifiers, annotations,
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
    map = Entity(
        entity_type=EntityTypeCv.PROTEIN,
        identifiers=Identifiers(
            # Primary UniProt accession
            Column('Entry', cv=IdentifierNamespaceCv.UNIPROT),

            # UniProt entry name
            Column('Entry Name', cv=IdentifierNamespaceCv.UNIPROT),

            # Primary gene name
            Column('Gene Names (primary)', cv=IdentifierNamespaceCv.GENE_NAME_PRIMARY),

            # Gene synonyms (space-delimited)
            Column(
                'Gene Names (synonym)',
                delimiter=' ',
                cv=IdentifierNamespaceCv.GENE_NAME_SYNONYM,
            ),

            # Protein name - extract primary (before parentheses)
            Column(
                'Protein names',
                processing={'extract_value': r'^([^(]+)'},
                cv=IdentifierNamespaceCv.NAME,
            ),

            # Protein name synonyms - split by '(' and extract content before ')'
            # Regex only matches tokens with ')', avoiding the primary name
            Column(
                'Protein names',
                delimiter='(',
                processing={'extract_value': r'^([^)]+)\)'},
                cv=IdentifierNamespaceCv.SYNONYM,
            ),

            # Cross-references (semicolon-delimited)
            Column('Ensembl', delimiter=';', cv=IdentifierNamespaceCv.ENSEMBL),
            Column('RefSeq', delimiter=';', cv=IdentifierNamespaceCv.REFSEQ),
            Column('PDB', delimiter=';', cv=IdentifierNamespaceCv.PDB),
            Column('AlphaFoldDB', delimiter=';', cv=IdentifierNamespaceCv.ALPHAFOLDDB),
            Column('KEGG', delimiter=';', cv=IdentifierNamespaceCv.KEGG),
            Column('ChEMBL', delimiter=';', cv=IdentifierNamespaceCv.CHEMBL),
            Column('SIGNOR', delimiter=';', cv=IdentifierNamespaceCv.SIGNOR),
            Column('IntAct', delimiter=';', cv=IdentifierNamespaceCv.INTACT),
            Column('BioGRID', delimiter=';', cv=IdentifierNamespaceCv.BIOGRID),
            Column('ComplexPortal', delimiter=';', cv=IdentifierNamespaceCv.COMPLEXPORTAL),
        ),
        annotations=Annotations(
            # Source annotation
            Column('Length', cv=MoleculeAnnotationsCv.SEQUENCE_LENGTH),
            Column('Mass', cv=MoleculeAnnotationsCv.MASS_DALTON),

            # Functional annotations
            Column('Function [CC]', cv=MoleculeAnnotationsCv.FUNCTION),
            Column('Subcellular location [CC]', cv=MoleculeAnnotationsCv.SUBCELLULAR_LOCATION),
            Column('Post-translational modification', cv=MoleculeAnnotationsCv.POST_TRANSLATIONAL_MODIFICATION),
            Column('Involvement in disease', cv=MoleculeAnnotationsCv.DISEASE_INVOLVEMENT),
            Column('Pathway', cv=MoleculeAnnotationsCv.PATHWAY_PARTICIPATION),
            Column('Activity regulation', cv=MoleculeAnnotationsCv.ACTIVITY_REGULATION),

            # Feature annotations
            Column('Mutagenesis', cv=MoleculeAnnotationsCv.MUTAGENESIS),
            Column('Transmembrane' ,cv=MoleculeAnnotationsCv.TRANSMEMBRANE_REGION),
            Column('Protein families', delimiter=",", cv=MoleculeAnnotationsCv.PROTEIN_FAMILY),
            Column('EC number', delimiter=";", cv=MoleculeAnnotationsCv.EC_NUMBER),

            # Sequence
            Column('Sequence', cv=MoleculeAnnotationsCv.AMINO_ACID_SEQUENCE),

            # GO terms as CV term accessions
            Column('Gene Ontology IDs', delimiter=';', cv=IdentifierNamespaceCv.CV_TERM_ACCESSION),

            # Keywords as CV term accessions
            Column('Keywords IDs', delimiter=';', cv=IdentifierNamespaceCv.CV_TERM_ACCESSION),

            # Organism as annotation
            Column('Organism (ID)', cv=IdentifierNamespaceCv.NCBI_TAX_ID),

            # PubMed references as annotations
            Column('PubMed ID', delimiter=';', cv=IdentifierNamespaceCv.PUBMED),
        ),
    )

    # Parse and yield entities
    if opener and opener.result:
        reader = csv.DictReader(opener.result, delimiter='\t')
        for row in reader:
            yield map(row)
