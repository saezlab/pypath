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
Parse Gene Ontology (GO) OBO data and emit Entity records.

This module converts Gene Ontology data into Entity records with CV_TERM type
using the schema defined in pypath.internals.silver_schema.
"""

from __future__ import annotations

from collections.abc import Generator

from pypath.formats.obo import Obo
from pypath.internals.silver_schema import Entity, Identifier, Annotation
from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    OntologyAnnotationCv,
    LicenseCV,
    UpdateCategoryCV,
    ResourceAnnotationCv,
    ResourceCv,
)
from ...internals.tabular_builder import (
    EntityBuilder,
    Annotations,
    Column,
    Identifiers,
)
from .shared import process_obo_term


# Gene Ontology OBO URL
GENE_ONTOLOGY_URL = "https://current.geneontology.org/ontology/go.obo"


def resource() -> Generator[Entity]:
    """
    Yield resource metadata as an Entity record.

    Yields:
        Entity record with type CV_TERM containing Gene Ontology metadata.
    """
    yield Entity(
        type=EntityTypeCv.CV_TERM,
        identifiers=[
            Identifier(type=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=ResourceCv.GENE_ONTOLOGY),
            Identifier(type=IdentifierNamespaceCv.NAME, value='Gene Ontology'),
        ],
        annotations=[
            Annotation(term=ResourceAnnotationCv.LICENSE, value=str(LicenseCV.CC_BY_4_0)),
            Annotation(term=ResourceAnnotationCv.UPDATE_CATEGORY, value=str(UpdateCategoryCV.REGULAR)),
            Annotation(term=IdentifierNamespaceCv.PUBMED, value='33290552'),
            Annotation(term=ResourceAnnotationCv.URL, value='http://geneontology.org/'),
            Annotation(term=ResourceAnnotationCv.DESCRIPTION, value=(
                'The Gene Ontology (GO) provides a comprehensive framework for '
                'describing gene functions across all organisms. GO covers three domains: '
                'biological process, molecular function, and cellular component. '
                'GO terms are organized in a hierarchical structure with rich relationships.'
            )),
        ],
    )


def gene_ontology() -> Generator[Entity]:
    """
    Download and parse Gene Ontology OBO file as Entity records.

    Downloads Gene Ontology OBO data and converts each term into a Entity
    with CV_TERM type, including identifiers, annotations, and relationships.

    Yields:
        Entity records with type CV_TERM
    """
    # Download and open the OBO file
    obo = Obo(GENE_ONTOLOGY_URL, name='GeneOntology')

    # Define the schema mapping
    schema = EntityBuilder(
        entity_type=EntityTypeCv.CV_TERM,
        identifiers=Identifiers(
            Column('accession', cv=IdentifierNamespaceCv.CV_TERM_ACCESSION),
            Column('name', cv=IdentifierNamespaceCv.NAME),
            Column('synonyms', delimiter=';', cv=IdentifierNamespaceCv.SYNONYM),
            Column('alt_ids', delimiter=';', cv=IdentifierNamespaceCv.CV_TERM_ACCESSION),

        ),
        annotations=Annotations(
            # Source annotation
            Column('definition', cv=OntologyAnnotationCv.DEFINITION),
            Column('comment', cv=OntologyAnnotationCv.COMMENT),
            Column('is_a', delimiter=';', cv=IdentifierNamespaceCv.CV_TERM_ACCESSION),
            Column('xrefs', delimiter=';', cv=IdentifierNamespaceCv.CV_TERM_ACCESSION),
            Column('is_obsolete', cv=OntologyAnnotationCv.IS_OBSOLETE),
        ),
    )

    # Parse and yield entities
    for term in obo:
        if term.stanza == 'Term':
            row = process_obo_term(term)
            yield schema(row)
