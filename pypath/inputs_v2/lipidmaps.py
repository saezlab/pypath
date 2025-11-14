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
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    Column,
    CV,
    EntityBuilder,
    IdentifiersBuilder,
    Map,
)
from pypath.share.downloads import download_and_open
from pypath.formats.sdf import SdfReader

def resource() -> Generator[Entity]:
    """
    Yield resource metadata as an Entity record.

    Yields:
        Entity record with type CV_TERM containing LIPID MAPS metadata.
    """
    yield Entity(
        type=EntityTypeCv.CV_TERM,
        identifiers=[
            Identifier(type=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=ResourceCv.LIPIDMAPS),
            Identifier(type=IdentifierNamespaceCv.NAME, value='LIPID MAPS Structure Database'),
        ],
        annotations=[
            Annotation(term=ResourceAnnotationCv.LICENSE, value=str(LicenseCV.CC_BY_4_0)),
            Annotation(term=ResourceAnnotationCv.UPDATE_CATEGORY, value=str(UpdateCategoryCV.REGULAR)),
            Annotation(term=IdentifierNamespaceCv.PUBMED, value='37855672'),
            Annotation(term=ResourceAnnotationCv.URL, value='https://lipidmaps.org/'),
            Annotation(term=ResourceAnnotationCv.DESCRIPTION, value=(
                'The LIPID MAPS Structure Database (LMSD) is a comprehensive database '
                'of lipid structures, annotations, and cross-references. It contains '
                'over 47,000 unique lipid structures, classified according to a '
                'comprehensive lipid classification system. The database includes '
                'structure-based identifiers (InChI, InChIKey, SMILES), chemical '
                'properties (formula, exact mass), lipid classification (category, '
                'main class, sub class), and cross-references to other databases '
                '(ChEBI, PubChem, SwissLipids, HMDB).'
            )),
        ],
    )


def lipidmaps_lipids(
    force_refresh: bool = False,
) -> Generator[Entity, None, None]:
    """
    Download and parse LIPID MAPS Structure Database as Entity records.

    This function downloads the LMSD structures in SDF format from the LIPID MAPS
    website, extracts the zip file, parses the SDF records, and yields Entity
    objects representing individual lipid structures.

    Args:
        force_refresh: If True, force redownload of the data.

    Yields:
        Entity records with type SMALL_MOLECULE, representing lipid structures
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

    def normalize_chebi(value: str | None) -> str | None:
        if not value:
            return None
        value = value.strip()
        if not value:
            return None
        return value if value.startswith('CHEBI:') else f'CHEBI:{value}'

    # Define declarative schema for LIPID MAPS lipids
    schema = EntityBuilder(
        entity_type=EntityTypeCv.LIPID,
        identifiers=IdentifiersBuilder(
            # Primary LIPID MAPS ID (LM_ID)
            CV(term=IdentifierNamespaceCv.LIPIDMAPS, value=Column('LM_ID')),
            # Common name
            CV(term=IdentifierNamespaceCv.NAME, value=Column('COMMON_NAME')),
            # Systematic name (IUPAC)
            CV(term=IdentifierNamespaceCv.IUPAC_NAME, value=Column('SYSTEMATIC_NAME')),
            # Abbreviation
            CV(term=IdentifierNamespaceCv.SYNONYM, value=Column('ABBREVIATION')),
            # Synonyms (may be semicolon-delimited)
            CV(term=IdentifierNamespaceCv.SYNONYM, value=Column('SYNONYMS', delimiter=';')),
            # Chemical structure identifiers
            CV(term=IdentifierNamespaceCv.STANDARD_INCHI_KEY, value=Column('INCHI_KEY')),
            CV(term=IdentifierNamespaceCv.STANDARD_INCHI, value=Column('INCHI')),
            CV(term=IdentifierNamespaceCv.SMILES, value=Column('SMILES')),
            # Chemical formula
            CV(term=IdentifierNamespaceCv.NAME, value=Column('FORMULA')),
            # Cross-references to other databases
            CV(
                term=IdentifierNamespaceCv.CHEBI,
                value=Map(col=Column('CHEBI_ID'), extract=[normalize_chebi]),
            ),
            CV(term=IdentifierNamespaceCv.PUBCHEM_COMPOUND, value=Column('PUBCHEM_CID')),
            CV(term=IdentifierNamespaceCv.HMDB, value=Column('HMDB_ID')),
            CV(term=IdentifierNamespaceCv.SWISSLIPIDS, value=Column('SWISSLIPIDS_ID')),
        ),
        annotations=AnnotationsBuilder(
            # Source annotation
            # Exact mass
            CV(term=MoleculeAnnotationsCv.MASS_DALTON, value=Column('EXACT_MASS')),
            # Lipid classification
            CV(term=MoleculeAnnotationsCv.LIPID_CATEGORY, value=Column('CATEGORY')),
            CV(term=MoleculeAnnotationsCv.LIPID_MAIN_CLASS, value=Column('MAIN_CLASS')),
            CV(term=MoleculeAnnotationsCv.LIPID_SUB_CLASS, value=Column('SUB_CLASS')),
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
