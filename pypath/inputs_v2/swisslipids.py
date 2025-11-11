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
Parse SwissLipids data and emit Entity records.

This module converts SwissLipids lipid data into Entity records using
the schema defined in pypath.internals.silver_schema.
"""

from __future__ import annotations

from collections.abc import Generator
import csv

from pypath.internals.silver_schema import Entity as SilverEntity, Resource
from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    LicenseCV,
    UpdateCategoryCV,
    MoleculeAnnotationsCv,
)
from pypath.share.downloads import download_and_open
from ..internals.tabular_builder import (
    Annotations,
    Column,
    Entity,
    Identifiers,
)


def get_resource() -> Resource:
    """
    Define the resource metadata.

    Returns:
        Resource object containing SwissLipids metadata.
    """
    return Resource(
        id='swisslipids',
        name='SwissLipids',
        license=LicenseCV.CC_BY_4_0,
        update_category=UpdateCategoryCV.REGULAR,
        publication='PMID:25943471',  # SwissLipids knowledgebase for lipid biology
        url='https://www.swisslipids.org/',
        description=(
            'SwissLipids is a curated resource providing a framework for the '
            'annotation of mass spectrometry data. It provides over 750,000 lipid '
            'structures with expert curation of lipid classes and nomenclature, '
            'hierarchical organization, cross-references to other databases (ChEBI, '
            'LIPID MAPS, HMDB), and integration with mass spectrometry tools. '
            'The database covers all major lipid categories including fatty acyls, '
            'glycerolipids, glycerophospholipids, sphingolipids, sterol lipids, '
            'prenol lipids, saccharolipids, and polyketides.'
        ),
    )


def swisslipids_lipids() -> Generator[SilverEntity, None, None]:
    """
    Download and parse SwissLipids lipid data as Entity records.

    Yields:
        Entity records with type LIPID, representing individual lipid molecules.
    """

    # Download the lipids table
    url = 'https://swisslipids.org/api/file.php?cas=download_files&file=lipids.tsv'
    opener = download_and_open(
        url=url,
        filename='lipids.tsv.gz',
        subfolder='swisslipids',
        encoding='latin-1',
    )

    # Define processing functions
    def filter_inchi(value):
        """Filter out placeholder InChI values."""
        return value and value != 'InChI=none'

    def filter_inchikey(value):
        """Filter out placeholder InChIKey values."""
        return value and value != 'InChIKey=none'

    def filter_chebi(value):
        """Format ChEBI IDs with CHEBI: prefix."""
        if not value or not value.strip():
            return None
        return f'CHEBI:{value}' if not value.startswith('CHEBI:') else value

    # Define the schema mapping
    schema = Entity(
        entity_type=EntityTypeCv.LIPID,
        identifiers=Identifiers(
            # Primary SwissLipids identifier
            Column('Lipid ID', cv=IdentifierNamespaceCv.SWISSLIPIDS),
            # Common name
            Column('Name', cv=IdentifierNamespaceCv.NAME),
            # Chemical structure identifiers
            Column('InChI key (pH7.3)', cv=IdentifierNamespaceCv.STANDARD_INCHI_KEY, processing={'filter': filter_inchikey}),
            Column('InChI (pH7.3)', cv=IdentifierNamespaceCv.STANDARD_INCHI, processing={'filter': filter_inchi}),
            Column('SMILES (pH7.3)', cv=IdentifierNamespaceCv.SMILES),
            # Cross-references to other databases
            Column('CHEBI', cv=IdentifierNamespaceCv.CHEBI, processing={'filter': filter_chebi}),
            Column('LIPID MAPS', cv=IdentifierNamespaceCv.LIPIDMAPS),
            Column('HMDB', cv=IdentifierNamespaceCv.HMDB),
            Column('MetaNetX', cv=IdentifierNamespaceCv.METANETX),
            # Synonyms (semicolon-separated)
            Column('Synonyms*', cv=IdentifierNamespaceCv.SYNONYM, delimiter=';'),
        ),
        annotations=Annotations(
            # Hierarchical classification
            Column('Level', cv=IdentifierNamespaceCv.CV_TERM_ACCESSION),
            Column('Lipid class*', cv=IdentifierNamespaceCv.NAME),
            Column('Parent', cv=IdentifierNamespaceCv.SWISSLIPIDS),
            Column('Components*', cv=IdentifierNamespaceCv.NAME),
            # Chemical properties
            Column('Formula (pH7.3)', cv=MoleculeAnnotationsCv.DESCRIPTION),
            Column('Charge (pH7.3)', cv=MoleculeAnnotationsCv.DESCRIPTION),
            Column('Exact Mass (neutral form)', cv=MoleculeAnnotationsCv.MASS_DALTON),
            # Display name
            Column('Abbreviation*', cv=IdentifierNamespaceCv.NAME),
            Column('PMID', delimiter='|', cv=IdentifierNamespaceCv.PUBMED)
        ),
    )

    # Parse and yield entities
    reader = csv.DictReader(opener.result, delimiter='\t')
    for row in reader:
        entity = schema(row)
        yield entity
