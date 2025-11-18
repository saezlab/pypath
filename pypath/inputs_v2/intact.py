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
Parse IntAct data and emit Entity records.

This module converts IntAct MITAB data into Entity records using the schema
defined in pypath.internals.silver_schema.
"""

from __future__ import annotations

from collections.abc import Generator

from pypath.share.downloads import download_and_open
from pypath.internals.silver_schema import Entity, Identifier, Annotation
from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    InteractionMetadataCv,
    ParticipantMetadataCv,
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
    Member,
    MembershipBuilder,
)
import csv


def resource() -> Generator[Entity]:
    """
    Yield resource metadata as an Entity record.

    Yields:
        Entity record with type CV_TERM containing IntAct metadata.
    """
    yield Entity(
        type=EntityTypeCv.CV_TERM,
        identifiers=[
            Identifier(type=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=ResourceCv.INTACT),
            Identifier(type=IdentifierNamespaceCv.NAME, value='IntAct'),
        ],
        annotations=[
            Annotation(term=ResourceAnnotationCv.LICENSE, value=str(LicenseCV.CC_BY_4_0)),
            Annotation(term=ResourceAnnotationCv.UPDATE_CATEGORY, value=str(UpdateCategoryCV.REGULAR)),
            Annotation(term=IdentifierNamespaceCv.PUBMED, value='37953288'),
            Annotation(term=ResourceAnnotationCv.URL, value='https://www.ebi.ac.uk/intact/'),
            Annotation(term=ResourceAnnotationCv.DESCRIPTION, value=(
                'IntAct provides a freely available, open source database system '
                'and analysis tools for molecular interaction data. All interactions '
                'are derived from literature curation or direct user submissions and '
                'are freely available in PSI-MITAB format. The database includes '
                'protein-protein, protein-small molecule and protein-nucleic acid '
                'interactions with detailed experimental evidence.'
            )),
        ],
    )


def intact_interactions(organism: int = 9606) -> Generator[Entity, None, None]:
    """
    Download and parse IntAct interactions as Entity records.

    Args:
        organism: NCBI taxonomy ID (9606 for human)

    Yields:
        Entity records with type INTERACTION, containing interactor pairs
    """
    if organism != 9606:
        raise ValueError(f'Currently only human (9606) is supported for IntAct')

    # Download and open the file
    url = 'https://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/species/human.zip'
    opener = download_and_open(
        url=url,
        filename='human.zip',
        subfolder='intact',
        large=True,
        ext='zip',
    )
    # Mapping of identifier prefixes (from IntAct MITAB data) to CV terms
    # All prefixes found in ID and Alt ID columns (from complete analysis of all 1.18M rows):
    # bind smid, cas registry number, chebi, chembl, chembl compound, ddbj/embl/genbank,
    # dip, ensembl, ensemblgenomes, entrezgene/locuslink, flybase, genbank identifier,
    # genbank_nucl_gi, genbank_protein_gi, hgnc, imex, intact, ipi, mint, mirbase,
    # pdbe, psi-mi, refseq, rfam, rnacentral, uniparc, uniprotkb
    identifier_cv_mapping = {
        'bind smid': IdentifierNamespaceCv.BIND,
        'cas registry number': IdentifierNamespaceCv.CAS,
        'chebi': IdentifierNamespaceCv.CHEBI,
        'chembl': IdentifierNamespaceCv.CHEMBL,
        'chembl compound': IdentifierNamespaceCv.CHEMBL_COMPOUND,
        'ddbj/embl/genbank': IdentifierNamespaceCv.REFSEQ,  # GenBank/EMBL/DDBJ accessions
        'dip': IdentifierNamespaceCv.DIP,
        'ensembl': IdentifierNamespaceCv.ENSEMBL,
        'ensemblgenomes': IdentifierNamespaceCv.ENSEMBL_GENOMES,
        'entrezgene/locuslink': IdentifierNamespaceCv.ENTREZ,
        'flybase': IdentifierNamespaceCv.FLYBASE,
        'genbank identifier': IdentifierNamespaceCv.GENBANK_IDENTIFIER,
        'genbank_nucl_gi': IdentifierNamespaceCv.GENBANK_NUCL_GI,
        'genbank_protein_gi': IdentifierNamespaceCv.GENBANK_PROTEIN_GI,
        'hgnc': IdentifierNamespaceCv.HGNC,
        'imex': IdentifierNamespaceCv.IMEX,
        'intact': IdentifierNamespaceCv.INTACT,
        'ipi': IdentifierNamespaceCv.IPI,
        'mint': IdentifierNamespaceCv.MINT,
        'mirbase': IdentifierNamespaceCv.MIRBASE,
        'pdbe': IdentifierNamespaceCv.PDB,  # PDBe is the European branch of PDB
        'psi-mi': IdentifierNamespaceCv.CV_TERM_ACCESSION,  # CV term accessions (e.g., MINT-xxx)
        'refseq': IdentifierNamespaceCv.REFSEQ,
        'rfam': IdentifierNamespaceCv.RFAM,
        'rnacentral': IdentifierNamespaceCv.RNACENTRAL,
        'uniparc': IdentifierNamespaceCv.UNIPARC,
        'uniprotkb': IdentifierNamespaceCv.UNIPROT,
    }

    prefix_regex = r'^([^:]+):'
    value_regex = r'^[^:]+:([^|"]+)'
    tax_regex = r'taxid:([-\d]+)'
    pubmed_regex = r'(?i)pubmed:(\d+)'
    mi_regex = r'(MI:\d+)'
    intact_value_regex = r'intact:([^|"]+)'

    def general_identifier_cv(column_name: str) -> CV:
        column = Column(column_name, delimiter='|')
        return CV(
            term=Map(
                col=column,
                extract=[prefix_regex, str.lower],
                map=identifier_cv_mapping,
            ),
            value=Map(col=column, extract=[value_regex]),
        )

    def mi_term_cv(column_name: str) -> CV:
        column = Column(column_name, delimiter='|')
        return CV(term=Map(col=column, extract=[mi_regex]))

    def tax_cv(column_name: str) -> CV:
        column = Column(column_name, delimiter='|')
        return CV(
            term=IdentifierNamespaceCv.NCBI_TAX_ID,
            value=Map(col=column, extract=[tax_regex]),
        )

    def pubmed_annotation(column_name: str) -> CV:
        column = Column(column_name, delimiter='|')
        return CV(
            term=IdentifierNamespaceCv.PUBMED,
            value=Map(col=column, extract=[pubmed_regex]),
        )

    def simple_annotation(term_cv, column_name: str) -> CV:
        return CV(term=term_cv, value=Column(column_name, delimiter='|'))

    # Define the schema mapping
    schema = EntityBuilder(
        entity_type=EntityTypeCv.INTERACTION,
        identifiers=IdentifiersBuilder(
            CV(
                term=IdentifierNamespaceCv.INTACT,
                value=Map(
                    col=Column('Interaction identifier(s)', delimiter='|'),
                    extract=[intact_value_regex],
                ),
            ),
        ),
        annotations=AnnotationsBuilder(
            # Source annotation
            # Interaction metadata
            mi_term_cv('Interaction type(s)'),
            mi_term_cv('Interaction detection method(s)'),
            mi_term_cv('Source database(s)'),
            simple_annotation(InteractionMetadataCv.CONFIDENCE_VALUE, 'Confidence value(s)'),
            simple_annotation(InteractionMetadataCv.EXPANSION_METHOD, 'Expansion method(s)'),

            # Publication information
            pubmed_annotation('Publication Identifier(s)'),
            # Experimental details
            tax_cv('Host organism(s)'),
            simple_annotation(InteractionMetadataCv.INTERACTION_PARAMETER, 'Interaction parameter(s)'),
            # Cross-references and annotations
            simple_annotation(InteractionMetadataCv.INTERACTION_XREF, 'Interaction Xref(s)'),
            simple_annotation(InteractionMetadataCv.INTERACTION_ANNOTATION, 'Interaction annotation(s)'),
            simple_annotation(InteractionMetadataCv.INTERACTION_CHECKSUM, 'Interaction Checksum(s)'),
        ),
        membership=MembershipBuilder(
            # Interactor A
            Member(
                entity=EntityBuilder(
                    entity_type=EntityTypeCv.PROTEIN,
                    identifiers=IdentifiersBuilder(
                        general_identifier_cv('#ID(s) interactor A'),
                        general_identifier_cv('Alt. ID(s) interactor A'),
                    ),
                    annotations=AnnotationsBuilder(
                        tax_cv('Taxid interactor A'),
                        simple_annotation(ParticipantMetadataCv.ALIAS, 'Alias(es) interactor A'),
                        simple_annotation(ParticipantMetadataCv.PARTICIPANT_XREF, 'Xref(s) interactor A'),
                        simple_annotation(ParticipantMetadataCv.PARTICIPANT_ANNOTATION, 'Annotation(s) interactor A'),
                        simple_annotation(ParticipantMetadataCv.PARTICIPANT_CHECKSUM, 'Checksum(s) interactor A'),
                    ),
                ),
                annotations=AnnotationsBuilder(
                    mi_term_cv('Biological role(s) interactor A'),
                    mi_term_cv('Experimental role(s) interactor A'),
                    mi_term_cv('Type(s) interactor A'),
                    simple_annotation(ParticipantMetadataCv.PARTICIPANT_FEATURE, 'Feature(s) interactor A'),
                    simple_annotation(ParticipantMetadataCv.STOICHIOMETRY, 'Stoichiometry(s) interactor A'),
                    mi_term_cv('Identification method participant A'),
                ),
            ),
            # Interactor B
            Member(
                entity=EntityBuilder(
                    entity_type=EntityTypeCv.PROTEIN,
                    identifiers=IdentifiersBuilder(
                        general_identifier_cv('ID(s) interactor B'),
                        general_identifier_cv('Alt. ID(s) interactor B'),
                    ),
                    annotations=AnnotationsBuilder(
                        tax_cv('Taxid interactor B'),
                        simple_annotation(ParticipantMetadataCv.ALIAS, 'Alias(es) interactor B'),
                        simple_annotation(ParticipantMetadataCv.PARTICIPANT_XREF, 'Xref(s) interactor B'),
                        simple_annotation(ParticipantMetadataCv.PARTICIPANT_ANNOTATION, 'Annotation(s) interactor B'),
                        simple_annotation(ParticipantMetadataCv.PARTICIPANT_CHECKSUM, 'Checksum(s) interactor B'),
                    ),
                ),
                annotations=AnnotationsBuilder(
                    mi_term_cv('Biological role(s) interactor B'),
                    mi_term_cv('Experimental role(s) interactor B'),
                    mi_term_cv('Type(s) interactor B'),
                    simple_annotation(ParticipantMetadataCv.PARTICIPANT_FEATURE, 'Feature(s) interactor B'),
                    simple_annotation(ParticipantMetadataCv.STOICHIOMETRY, 'Stoichiometry(s) interactor B'),
                    mi_term_cv('Identification method participant B'),
                ),
            ),
        ),
    )

    # Get the file from the zip (it's a dict with filename -> file handle)
    for file_handle in opener.result.values():
        # Parse and yield entities
        reader = csv.DictReader(file_handle, delimiter='\t')
        for row in reader:
            yield schema(row)
        break  # Only process the first file in the zip
