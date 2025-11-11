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
Parse UniProt Keywords OBO data and emit Entity records.

This module converts UniProt Keywords ontology data into Entity records with CV_TERM type
using the schema defined in pypath.internals.silver_schema.
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


# UniProt Keywords OBO URL
UNIPROT_KEYWORDS_URL = "https://rest.uniprot.org/keywords/stream?format=obo&query=%28*%29"


def get_resource() -> Resource:
    """
    Define the resource metadata.

    Returns:
        Resource object containing UniProt Keywords metadata.
    """
    return Resource(
        id='uniprot_keywords',
        name='UniProt Keywords',
        license=LicenseCV.CC_BY_4_0,
        update_category=UpdateCategoryCV.REGULAR,
        publication='PMID:33237286',
        url='https://www.uniprot.org/keywords/',
        description=(
            'UniProt Keywords provide a controlled vocabulary for '
            'summarizing protein properties, including biological processes, '
            'molecular functions, cellular components, protein families, and more. '
            'Keywords facilitate consistent annotation and efficient searching.'
        ),
    )


def uniprot_keywords() -> Generator[SilverEntity]:
    """
    Download and parse UniProt Keywords OBO file as Entity records.

    Downloads UniProt Keywords OBO data and converts each term into a SilverEntity
    with CV_TERM type, including identifiers, annotations, and relationships.

    Yields:
        Entity records with type CV_TERM
    """
    # Download and open the OBO file
    obo = Obo(UNIPROT_KEYWORDS_URL, name='UniProtKeywords')

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
            Column('comment', cv=OntologyAnnotationCv.COMMENT),
            Column('is_a', delimiter=';', cv=IdentifierNamespaceCv.CV_TERM_ACCESSION),
            Column('xrefs', delimiter=';', cv=IdentifierNamespaceCv.CV_TERM_ACCESSION),

        ),
    )

    # Parse and yield entities
    for term in obo:
        if term.stanza == 'Term':
            row = process_obo_term(term)
            yield schema(row)
