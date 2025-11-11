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
Parse Guide to Pharmacology data and emit Entity records.

This module converts Guide to Pharmacology ligand-target interaction data
into Entity records using the schema defined in pypath.internals.silver_schema.
"""

from __future__ import annotations

from collections.abc import Generator
import csv

from pypath.share.downloads import download_and_open
from pypath.internals.silver_schema import Entity as SilverEntity, Resource
from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    LicenseCV,
    UpdateCategoryCV,
    InteractionParameterCv,
    AnnotationTypeCv,
)
from ..internals.tabular_builder import (
    Annotations,
    Column,
    Entity,
    Identifiers,
    Member,
    Members,
)


def get_resource() -> Resource:
    """
    Define the resource metadata.

    Returns:
        Resource object containing Guide to Pharmacology metadata.
    """
    return Resource(
        id='guidetopharma',
        name='Guide to Pharmacology',
        license=LicenseCV.CC_BY_SA_4_0,
        update_category=UpdateCategoryCV.REGULAR,
        publication='PMID:37953350',  # GtoPdb in 2024 paper
        url='https://www.guidetopharmacology.org/',
        description=(
            'The IUPHAR/BPS Guide to PHARMACOLOGY is an expert-curated resource '
            'of ligand-activity-target relationships, providing quantitative '
            'information on drug targets and the prescription medicines and '
            'experimental drugs that act on them. It covers G protein-coupled '
            'receptors, voltage-gated ion channels, ligand-gated ion channels, '
            'nuclear hormone receptors, catalytic receptors, enzymes, and '
            'transporters.'
        ),
    )


def guidetopharma_interactions() -> Generator[SilverEntity, None, None]:
    """
    Download and parse Guide to Pharmacology ligand-target interactions as Entity records.

    Yields:
        Entity records with type INTERACTION, representing ligand-target interactions.
    """

    # Download the interactions table
    url = 'http://www.guidetopharmacology.org/DATA/interactions.csv'
    opener = download_and_open(
        url=url,
        filename='interactions.csv',
        subfolder='guidetopharma',
    )

    # Define the schema mapping
    schema = Entity(
        entity_type=EntityTypeCv.INTERACTION,
        identifiers=Identifiers(
            # Interaction identifiers
            Column('Interaction ID', cv=IdentifierNamespaceCv.GUIDETOPHARMA),
        ),
        annotations=Annotations(
            # Interaction properties
            Column('Action', cv=AnnotationTypeCv.ACTION),
            Column('Type', cv=AnnotationTypeCv.ACTION_TYPE),
            Column('Endogenous', cv=AnnotationTypeCv.ENDOGENOUS),
            Column('Affinity High', cv=InteractionParameterCv.AFFINITY_HIGH),
            Column('Affinity Low', cv=InteractionParameterCv.AFFINITY_LOW),
            Column('Affinity Median', cv=InteractionParameterCv.AFFINITY_MEDIAN),
            Column('Affinity Units', cv=AnnotationTypeCv.AFFINITY_UNITS),
            Column('Primary Target', cv=AnnotationTypeCv.PRIMARY_TARGET),
            Column('PubMed ID', cv=IdentifierNamespaceCv.PUBMED),
            Column('Target Species', cv=AnnotationTypeCv.ORGANISM),
        ),
        members=Members(
            # Ligand
            Member(
                entity=Entity(
                    entity_type=EntityTypeCv.SMALL_MOLECULE,
                    identifiers=Identifiers(
                        Column('Ligand ID', cv=IdentifierNamespaceCv.GUIDETOPHARMA),
                    ),
                ),
            ),
            # Target
            Member(
                entity=Entity(
                    entity_type=EntityTypeCv.PROTEIN,
                    identifiers=Identifiers(
                        Column('Target ID', cv=IdentifierNamespaceCv.GUIDETOPHARMA),
                    ),
                ),
            ),
        ),
    )

    # Parse and yield entities
    # Skip the first line (version info)
    next(opener.result)
    reader = csv.DictReader(opener.result)
    for row in reader:
        yield schema(row)
