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
Parse PSI-MI Controlled Vocabulary OBO data and emit Entity records.

This module converts PSI-MI (Proteomics Standards Initiative - Molecular Interactions)
controlled vocabulary data into Entity records with CV_TERM type using the schema
defined in pypath.internals.silver_schema.
"""

from __future__ import annotations

from collections.abc import Generator

from pypath.formats.obo import Obo
from pypath.internals.silver_schema import Entity as SilverEntity, Resource
from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    OntologyAnnotationCv,
    LicenseCV,
    UpdateCategoryCV,
)
from ...internals.tabular_builder import (
    Annotations,
    Column,
    Entity,
    Identifiers,
)
from .shared import process_obo_term


# PSI-MI OBO URL
PSI_MI_URL = "https://raw.githubusercontent.com/HUPO-PSI/psi-mi-CV/master/psi-mi.obo"


def get_resource() -> Resource:
    """
    Define the resource metadata.

    Returns:
        Resource object containing PSI-MI metadata.
    """
    return Resource(
        id='psi_mi',
        name='PSI-MI Controlled Vocabulary',
        license=LicenseCV.CC_BY_4_0,
        update_category=UpdateCategoryCV.REGULAR,
        publication='PMID:17925023',  # PSI-MI 2.5 paper
        url='https://github.com/HUPO-PSI/psi-mi-CV',
        description=(
            'The PSI-MI (Proteomics Standards Initiative - Molecular Interactions) '
            'controlled vocabulary provides standardized terms for describing '
            'molecular interactions. It includes terms for interaction types, '
            'detection methods, participant roles, and experimental features.'
        ),
    )


def psi_mi() -> Generator[SilverEntity]:
    """
    Download and parse PSI-MI OBO file as Entity records.

    Downloads PSI-MI OBO data and converts each term into a SilverEntity
    with CV_TERM type, including identifiers, annotations, and relationships.

    Yields:
        Entity records with type CV_TERM
    """
    # Download and open the OBO file
    obo = Obo(PSI_MI_URL, name='PsiMi')

    # Define the schema mapping
    schema = Entity(
        entity_type=EntityTypeCv.CV_TERM,
        identifiers=Identifiers(
            Column('accession', cv=IdentifierNamespaceCv.CV_TERM_ACCESSION),
            Column('name', cv=IdentifierNamespaceCv.NAME),
            Column('synonyms', delimiter=';', cv=IdentifierNamespaceCv.SYNONYM),
            Column('alt_ids', delimiter=';', cv=IdentifierNamespaceCv.CV_TERM_ACCESSION),

        ),
        annotations=Annotations(
            Column('definition', cv=OntologyAnnotationCv.DEFINITION),
            Column('is_a', delimiter=';', cv=IdentifierNamespaceCv.CV_TERM_ACCESSION),
            Column('xrefs', delimiter=';', cv=IdentifierNamespaceCv.CV_TERM_ACCESSION),
            Column('comment', cv=OntologyAnnotationCv.COMMENT),
            Column('is_obsolete', cv=OntologyAnnotationCv.IS_OBSOLETE),
        ),
    )

    # Parse and yield entities
    for term in obo:
        if term.stanza == 'Term':
            row = process_obo_term(term)
            yield schema(row)
