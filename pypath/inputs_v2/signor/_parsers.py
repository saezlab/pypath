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
Parse SIGNOR data and emit Entity records.

This module converts SIGNOR data into Entity records using the schema
defined in pypath.internals.silver_schema.
"""

from __future__ import annotations

from collections.abc import Generator

from pypath.share.downloads import download_and_open
from pypath.internals.silver_schema import Entity as SilverEntity
from pypath.internals.cv_terms import EntityTypeCv, IdentifierNamespaceCv
from ..tabular_builder import (
    Annotations,
    Column,
    Entity as EntitySchema,
    Entities as MemberEntities,
    Identifiers,
    Members,
)
import csv
from pypath.internals.silver_schema import Resource
from omnipath_build.utils.cv_terms import LicenseCV, UpdateCategoryCV

__all__ = ['SIGNOR_RESOURCE']


def get_resource() -> Resource:
    """
    Define the resource metadata.

    Returns:
        Resource object containing SIGNOR metadata.
    """
    return Resource(
        id='signor',
        name='SIGNOR',
        license=LicenseCV.CC_BY_4_0,
        update_category=UpdateCategoryCV.REGULAR,
        publication='PMID:31665520',  # SIGNOR 2.0 paper
        url='https://signor.uniroma2.it/',
        description=(
            'SIGNOR (SIGnaling Network Open Resource) is a comprehensive '
            'resource of causal relationships between biological entities '
            'with a focus on signaling pathways. It provides manually curated '
            'interactions with mechanistic details including protein-protein '
            'interactions, post-translational modifications, transcriptional '
            'regulation, and small molecule effects.'
        ),
    )


def signor_complexes() -> Generator[SilverEntity]:
    """
    Download and parse SIGNOR complex data as Entity records.

    Yields:
        Entity records with type PROTEIN_COMPLEX, containing member proteins
    """
    # Download and open the file
    url = 'https://signor.uniroma2.it/download_complexes.php'
    opener = download_and_open(
        url,
        filename='signor_complexes.txt',
        subfolder='signor',
        query={'submit': 'Download complex data'},
        post=True,
    )

    # Define the schema mapping
    map = EntitySchema(
        entity_type=EntityTypeCv.PROTEIN_COMPLEX,
        identifiers=Identifiers(
            Column('SIGNOR ID', cv=IdentifierNamespaceCv.SIGNOR),
            Column('COMPLEX NAME', cv=IdentifierNamespaceCv.NAME),
        ),
        members=Members(
            MemberEntities(
                entity_type=EntityTypeCv.PROTEIN,
                identifiers=Identifiers(
                    Column('LIST OF ENTITIES', delimiter=',', cv=IdentifierNamespaceCv.UNIPROT),
                ),
            )
        ),
    )

    # Parse and yield entities
    for entity in csv.DictReader(opener.result, delimiter=';'):
        yield map(entity)


def signor_protein_families() -> Generator[SilverEntity]:
    """
    Download and parse SIGNOR protein family data as Entity records.

    Yields:
        Entity records with type PROTEIN_FAMILY, containing member proteins
    """
    # Download and open the file
    url = 'https://signor.uniroma2.it/download_complexes.php'
    opener = download_and_open(
        url,
        filename='signor_protein_families.txt',
        subfolder='signor',
        query={'submit': 'Download protein family data'},
        post=True,
    )

    # Define the schema mapping
    map = EntitySchema(
        entity_type=EntityTypeCv.PROTEIN_FAMILY,
        identifiers=Identifiers(
            Column('SIGNOR ID', cv=IdentifierNamespaceCv.SIGNOR),
            Column('PROT. FAMILY NAME', cv=IdentifierNamespaceCv.NAME),
        ),
        members=Members(
            MemberEntities(
                entity_type=EntityTypeCv.PROTEIN,
                identifiers=Identifiers(
                    Column('LIST OF ENTITIES', delimiter=',', cv=IdentifierNamespaceCv.UNIPROT),
                ),
            )
        ),
    )

    # Parse and yield entities
    for entity in csv.DictReader(opener.result, delimiter=';'):
        yield map(entity)


def signor_phenotypes() -> Generator[SilverEntity]:
    """
    Download and parse SIGNOR phenotype data as Entity records.

    Yields:
        Entity records with type PHENOTYPE
    """
    # Download and open the file
    url = 'https://signor.uniroma2.it/download_complexes.php'
    opener = download_and_open(
        url,
        filename='signor_phenotypes.txt',
        subfolder='signor',
        query={'submit': 'Download phenotype data'},
        post=True,
    )

    # Define the schema mapping
    map = EntitySchema(
        entity_type=EntityTypeCv.PHENOTYPE,
        identifiers=Identifiers(
            Column('SIGNOR ID', cv=IdentifierNamespaceCv.SIGNOR),
            Column('PHENOTYPE NAME', cv=IdentifierNamespaceCv.NAME),
        ),
        annotations=Annotations(
            Column('PHENOTYPE DESCRIPTION', cv=IdentifierNamespaceCv.SYNONYM),
        ),
    )

    # Parse and yield entities
    for entity in csv.DictReader(opener.result, delimiter=';'):
        yield map(entity)


def signor_stimuli() -> Generator[SilverEntity]:
    """
    Download and parse SIGNOR stimulus data as Entity records.

    Yields:
        Entity records with type STIMULUS
    """
    # Download and open the file
    url = 'https://signor.uniroma2.it/download_complexes.php'
    opener = download_and_open(
        url,
        filename='signor_stimuli.txt',
        subfolder='signor',
        query={'submit': 'Download stimulus data'},
        post=True,
    )

    # Define the schema mapping
    map = EntitySchema(
        entity_type=EntityTypeCv.STIMULUS,
        identifiers=Identifiers(
            Column('SIGNOR ID', cv=IdentifierNamespaceCv.SIGNOR),
            Column('STIMULUS NAME', cv=IdentifierNamespaceCv.NAME),
        ),
        annotations=Annotations(
            Column('STIMULUS DESCRIPTION', cv=IdentifierNamespaceCv.SYNONYM),
        ),
    )

    # Parse and yield entities
    for entity in csv.DictReader(opener.result, delimiter=';'):
        yield map(entity)
