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
from pypath.internals.cv_terms import EntityTypeCv, IdentifierNamespaceCv, LicenseCV, UpdateCategoryCV, ResourceAnnotationCv, ResourceCv
from ..internals.tabular_builder import (
    EntityBuilder,
    Annotations,
    Column,
    Identifiers,
    Member,
    Members,
)
import csv


def intact() -> Generator[Entity]:
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

    # Processing patterns for MITAB fields
    general_id_processing = {'extract_prefix': r'^([^:]+):', 'extract_value': r'^[^:]+:([^|"]+)'}
    tax_processing = {'extract_value': r'taxid:([-\d]+)'}
    pubmed_processing = {'extract_value': r'(?i)pubmed:(\d+)'}
    mi_term_processing = {'extract_term': r'(MI:\d+)'}
    intact_processing = {'extract_value': r'intact:([^|"]+)'}

    # Define the schema mapping
    schema = EntityBuilder(
        entity_type=EntityTypeCv.INTERACTION,
        identifiers=Identifiers(
            Column(
                'Interaction identifier(s)',
                delimiter='|',
                processing=intact_processing,
                cv=IdentifierNamespaceCv.INTACT,
            ),
        ),
        annotations=Annotations(
            # Source annotation
            # Interaction metadata
            Column('Interaction type(s)', delimiter='|', processing=mi_term_processing),
            Column('Interaction detection method(s)', delimiter='|', processing=mi_term_processing),
            Column('Source database(s)', delimiter='|', processing=mi_term_processing),
            Column('Confidence value(s)', delimiter='|'),
            Column('Expansion method(s)', delimiter='|'),

            # Publication information
            Column('Publication Identifier(s)', delimiter='|', processing=pubmed_processing, cv=IdentifierNamespaceCv.PUBMED),
            Column('Publication 1st author(s)', delimiter='|'),

            # Experimental details
            Column('Host organism(s)', delimiter='|', processing=tax_processing),
            Column('Interaction parameter(s)', delimiter='|'),
            Column('Negative', delimiter='|'),

            # Cross-references and annotations
            Column('Interaction Xref(s)', delimiter='|'),
            Column('Interaction annotation(s)', delimiter='|'),
            Column('Interaction Checksum(s)', delimiter='|'),

            # Timestamps
            Column('Creation date', delimiter='|'),
            Column('Update date', delimiter='|'),
        ),
        members=Members(
            # Interactor A
            Member(
                entity=EntityBuilder(
                    entity_type=EntityTypeCv.PROTEIN,
                    identifiers=Identifiers(
                        Column('#ID(s) interactor A', delimiter='|', processing=general_id_processing, cv=identifier_cv_mapping),
                        Column('Alt. ID(s) interactor A', delimiter='|', processing=general_id_processing, cv=identifier_cv_mapping),
                    ),
                    annotations=Annotations(
                        Column('Taxid interactor A', delimiter='|', processing=tax_processing, cv=IdentifierNamespaceCv.NCBI_TAX_ID),
                        Column('Alias(es) interactor A', delimiter='|'),
                        Column('Xref(s) interactor A', delimiter='|'),
                        Column('Annotation(s) interactor A', delimiter='|'),
                        Column('Checksum(s) interactor A', delimiter='|'),
                    ),
                ),
                annotations=Annotations(
                    Column('Biological role(s) interactor A', delimiter='|', processing=mi_term_processing),
                    Column('Experimental role(s) interactor A', delimiter='|', processing=mi_term_processing),
                    Column('Type(s) interactor A', delimiter='|', processing=mi_term_processing),
                    Column('Feature(s) interactor A', delimiter='|'),
                    Column('Stoichiometry(s) interactor A', delimiter='|'),
                    Column('Identification method participant A', delimiter='|', processing=mi_term_processing),
                ),
            ),
            # Interactor B
            Member(
                entity=EntityBuilder(
                    entity_type=EntityTypeCv.PROTEIN,
                    identifiers=Identifiers(
                        Column('ID(s) interactor B', delimiter='|', processing=general_id_processing, cv=identifier_cv_mapping),
                        Column('Alt. ID(s) interactor B', delimiter='|', processing=general_id_processing, cv=identifier_cv_mapping),
                    ),
                    annotations=Annotations(
                        Column('Taxid interactor B', delimiter='|', processing=tax_processing, cv=IdentifierNamespaceCv.NCBI_TAX_ID),
                        Column('Alias(es) interactor B', delimiter='|'),
                        Column('Xref(s) interactor B', delimiter='|'),
                        Column('Annotation(s) interactor B', delimiter='|'),
                        Column('Checksum(s) interactor B', delimiter='|'),
                    ),
                ),
                annotations=Annotations(
                    Column('Biological role(s) interactor B', delimiter='|', processing=mi_term_processing),
                    Column('Experimental role(s) interactor B', delimiter='|', processing=mi_term_processing),
                    Column('Type(s) interactor B', delimiter='|', processing=mi_term_processing),
                    Column('Feature(s) interactor B', delimiter='|'),
                    Column('Stoichiometry(s) interactor B', delimiter='|'),
                    Column('Identification method participant B', delimiter='|', processing=mi_term_processing),
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
