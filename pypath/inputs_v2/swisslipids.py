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

from pypath.internals.silver_schema import Entity, Identifier, Annotation
from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    LicenseCV,
    UpdateCategoryCV,
    MoleculeAnnotationsCv,
    ResourceAnnotationCv,
    ResourceCv,
)
from pypath.share.downloads import download_and_open
from ..internals.tabular_builder import (
    AnnotationsBuilder,
    Column,
    CV,
    EntityBuilder,
    IdentifiersBuilder,
    Map,
)


def resource() -> Generator[Entity]:
    """
    Yield resource metadata as an Entity record.

    Yields:
        Entity record with type CV_TERM containing SwissLipids metadata.
    """
    yield Entity(
        type=EntityTypeCv.CV_TERM,
        identifiers=[
            Identifier(type=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=ResourceCv.SWISSLIPIDS),
            Identifier(type=IdentifierNamespaceCv.NAME, value='SwissLipids'),
        ],
        annotations=[
            Annotation(term=ResourceAnnotationCv.LICENSE, value=str(LicenseCV.CC_BY_4_0)),
            Annotation(term=ResourceAnnotationCv.UPDATE_CATEGORY, value=str(UpdateCategoryCV.REGULAR)),
            Annotation(term=IdentifierNamespaceCv.PUBMED, value='25943471'),
            Annotation(term=ResourceAnnotationCv.URL, value='https://www.swisslipids.org/'),
            Annotation(term=ResourceAnnotationCv.DESCRIPTION, value=(
                'SwissLipids is a curated resource providing a framework for the '
                'annotation of mass spectrometry data. It provides over 750,000 lipid '
                'structures with expert curation of lipid classes and nomenclature, '
                'hierarchical organization, cross-references to other databases (ChEBI, '
                'LIPID MAPS, HMDB), and integration with mass spectrometry tools. '
                'The database covers all major lipid categories including fatty acyls, '
                'glycerolipids, glycerophospholipids, sphingolipids, sterol lipids, '
                'prenol lipids, saccharolipids, and polyketides.'
            )),
        ],
    )


def swisslipids_lipids() -> Generator[Entity, None, None]:
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
        if not value or value == 'InChI=none':
            return None
        return value

    def filter_inchikey(value):
        """Filter out placeholder InChIKey values."""
        if not value or value == 'InChIKey=none':
            return None
        return value

    def filter_chebi(value):
        """Format ChEBI IDs with CHEBI: prefix."""
        if not value or not value.strip():
            return None
        return f'CHEBI:{value}' if not value.startswith('CHEBI:') else value

    # Define the schema mapping
    schema = EntityBuilder(
        entity_type=EntityTypeCv.LIPID,
        identifiers=IdentifiersBuilder(
            # Primary SwissLipids identifier
            CV(term=IdentifierNamespaceCv.SWISSLIPIDS, value=Column('Lipid ID')),
            CV(term=IdentifierNamespaceCv.NAME, value=Column('Name')),
            CV(
                term=IdentifierNamespaceCv.STANDARD_INCHI_KEY,
                value=Map(col=Column('InChI key (pH7.3)'), extract=[filter_inchikey]),
            ),
            CV(
                term=IdentifierNamespaceCv.STANDARD_INCHI,
                value=Map(col=Column('InChI (pH7.3)'), extract=[filter_inchi]),
            ),
            CV(term=IdentifierNamespaceCv.SMILES, value=Column('SMILES (pH7.3)')),
            CV(
                term=IdentifierNamespaceCv.CHEBI,
                value=Map(col=Column('CHEBI'), extract=[filter_chebi]),
            ),
            CV(term=IdentifierNamespaceCv.LIPIDMAPS, value=Column('LIPID MAPS')),
            CV(term=IdentifierNamespaceCv.HMDB, value=Column('HMDB')),
            CV(term=IdentifierNamespaceCv.METANETX, value=Column('MetaNetX')),
            CV(term=IdentifierNamespaceCv.SYNONYM, value=Column('Synonyms*', delimiter=';')),
            CV(term=IdentifierNamespaceCv.SYNONYM, value=Column('Abbreviation*')),
        ),
        annotations=AnnotationsBuilder(
            # Source annotation
            CV(term=MoleculeAnnotationsCv.LIPID_HIERARCHY_LEVEL, value=Column('Level')),
            CV(term=MoleculeAnnotationsCv.LIPID_MAIN_CLASS, value=Column('Lipid class*')),
            CV(term=IdentifierNamespaceCv.SWISSLIPIDS, value=Column('Parent')),
            CV(term=MoleculeAnnotationsCv.LIPID_STRUCTURAL_COMPONENTS, value=Column('Components*')),
            CV(term=MoleculeAnnotationsCv.MOLECULAR_CHARGE, value=Column('Charge (pH7.3)')),
            CV(term=MoleculeAnnotationsCv.MASS_DALTON, value=Column('Exact Mass (neutral form)')),
            CV(term=IdentifierNamespaceCv.PUBMED, value=Column('PMID', delimiter='|')),
        ),
    )

    # Parse and yield entities
    reader = csv.DictReader(opener.result, delimiter='\t')
    for row in reader:
        entity = schema(row)
        yield entity
