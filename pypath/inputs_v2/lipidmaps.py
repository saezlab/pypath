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
Parse LIPID MAPS Structure Database (LMSD) data and emit Entity records.

This module downloads the LIPID MAPS Structure Database in SDF format and
converts lipid structures into Entity records using the schema defined in
pypath.internals.silver_schema.
"""

from __future__ import annotations

from collections.abc import Generator

from pypath.internals.silver_schema import Entity as SilverEntity, Resource
from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    MoleculeAnnotationsCv,
    LicenseCV,
    UpdateCategoryCV,
)
from pypath.internals.tabular_builder import (
    Annotations,
    Column,
    Entity,
    Identifiers,
)
from pypath.share.downloads import download_and_open
from pypath.formats.sdf import SdfReader


def get_resource() -> Resource:
    """
    Define the resource metadata.

    Returns:
        Resource object containing LIPID MAPS metadata.
    """
    return Resource(
        id='lipidmaps',
        name='LIPID MAPS Structure Database',
        license=LicenseCV.CC_BY_4_0,
        update_category=UpdateCategoryCV.REGULAR,
        publication='PMID:37855672',  # LIPID MAPS 2024 update paper
        url='https://lipidmaps.org/',
        description=(
            'The LIPID MAPS Structure Database (LMSD) is a comprehensive database '
            'of lipid structures, annotations, and cross-references. It contains '
            'over 47,000 unique lipid structures, classified according to a '
            'comprehensive lipid classification system. The database includes '
            'structure-based identifiers (InChI, InChIKey, SMILES), chemical '
            'properties (formula, exact mass), lipid classification (category, '
            'main class, sub class), and cross-references to other databases '
            '(ChEBI, PubChem, SwissLipids, HMDB).'
        ),
    )


def lipidmaps_lipids(
    force_refresh: bool = False,
) -> Generator[SilverEntity, None, None]:
    """
    Download and parse LIPID MAPS Structure Database as Entity records.

    This function downloads the LMSD structures in SDF format from the LIPID MAPS
    website, extracts the zip file, parses the SDF records, and yields SilverEntity
    objects representing individual lipid structures.

    Args:
        force_refresh: If True, force redownload of the data.

    Yields:
        SilverEntity records with type SMALL_MOLECULE, representing lipid structures
        with their identifiers (LIPID MAPS ID, InChI, SMILES, etc.), chemical
        properties (formula, mass, classification), and cross-references to other
        databases (ChEBI, PubChem, SwissLipids, HMDB).
    """
    # Download LIPID MAPS structure database SDF file
    url = 'https://lipidmaps.org/files/?file=LMSD&ext=sdf.zip'
    opener = download_and_open(
        url,
        filename='structures.zip',
        subfolder='lipidmaps',
        large=True,
        ext='zip',
        default_mode='rb',  # Binary mode required for SDF parser
        force=force_refresh,
    )

    # Define declarative schema for LIPID MAPS lipids
    schema = Entity(
        entity_type=EntityTypeCv.SMALL_MOLECULE,
        identifiers=Identifiers(
            # Primary LIPID MAPS ID (LM_ID)
            Column('LM_ID', cv=IdentifierNamespaceCv.LIPIDMAPS),
            # Common name
            Column('COMMON_NAME', cv=IdentifierNamespaceCv.NAME),
            # Systematic name (IUPAC)
            Column('SYSTEMATIC_NAME', cv=IdentifierNamespaceCv.IUPAC_NAME),
            # Abbreviation
            Column('ABBREVIATION', cv=IdentifierNamespaceCv.SYNONYM),
            # Synonyms (may be semicolon-delimited)
            Column('SYNONYMS', delimiter=';', cv=IdentifierNamespaceCv.SYNONYM),
            # Chemical structure identifiers
            Column('INCHI_KEY', cv=IdentifierNamespaceCv.STANDARD_INCHI_KEY),
            Column('INCHI', cv=IdentifierNamespaceCv.STANDARD_INCHI),
            Column('SMILES', cv=IdentifierNamespaceCv.SMILES),
            # Chemical formula
            Column('FORMULA', cv=IdentifierNamespaceCv.NAME),
            # Cross-references to other databases
            Column(
                'CHEBI_ID',
                cv=IdentifierNamespaceCv.CHEBI,
                processing={
                    'extract_prefix': r'^(CHEBI:)?',
                    'extract_value': r'^(?:CHEBI:)?(.+)',
                },
            ),
            Column('PUBCHEM_CID', cv=IdentifierNamespaceCv.PUBCHEM_COMPOUND),
            Column('HMDB_ID', cv=IdentifierNamespaceCv.HMDB),
            Column('SWISSLIPIDS_ID', cv=IdentifierNamespaceCv.SWISSLIPIDS),
        ),
        annotations=Annotations(
            # Exact mass
            Column('EXACT_MASS', cv=MoleculeAnnotationsCv.MASS_DALTON),
            # Lipid classification
            Column('CATEGORY', cv=MoleculeAnnotationsCv.DESCRIPTION),
            Column('MAIN_CLASS', cv=MoleculeAnnotationsCv.DESCRIPTION),
            Column('SUB_CLASS', cv=MoleculeAnnotationsCv.DESCRIPTION),
            # Additional names and descriptions
            Column('NAME', cv=MoleculeAnnotationsCv.DESCRIPTION),
        ),
    )

    # Parse and yield entities
    if opener and opener.result:
        # Extract the SDF file from the zip
        for file_handle in opener.result.values():
            # Use SdfReader to parse the SDF file
            sdf_reader = SdfReader(
                file_handle,
                names={
                    'HMDB_ID': 'HMDB_ID',
                    'PUBCHEM_CID': 'PUBCHEM_CID',
                    'SWISSLIPIDS_ID': 'SWISSLIPIDS_ID',
                    'LM_ID': 'LM_ID',
                    'ABBREVIATION': 'ABBREVIATION',
                    'CHEBI_ID': 'CHEBI_ID',
                    'SYNONYMS': 'SYNONYMS',
                    'INCHI': 'INCHI',
                    'INCHI_KEY': 'INCHI_KEY',
                    'COMMON_NAME': 'COMMON_NAME',
                    'SYSTEMATIC_NAME': 'SYSTEMATIC_NAME',
                    'SMILES': 'SMILES',
                    'FORMULA': 'FORMULA',
                },
                fields={
                    'EXACT_MASS',
                    'CATEGORY',
                    'MAIN_CLASS',
                    'SUB_CLASS',
                    'NAME',
                },
            )

            for record in sdf_reader:
                # Flatten the record structure to match the schema
                # The SDF reader returns records with 'name' and 'annot' dicts
                flat_record = {}

                # Copy identifiers from 'name' dict
                if 'name' in record:
                    flat_record.update(record['name'])

                # Copy annotations from 'annot' dict
                if 'annot' in record:
                    flat_record.update(record['annot'])

                # Yield the entity
                yield schema(flat_record)

            # Only process the first file in the zip
            break
