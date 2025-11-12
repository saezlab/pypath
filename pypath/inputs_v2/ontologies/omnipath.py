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
Extract OmniPath controlled vocabulary terms and their parent relationships.

This module converts OmniPath CV terms from pypath.internals.cv_terms into Entity records
with CV_TERM type, including their hierarchical parent relationships (is_a).
"""

from __future__ import annotations

from collections.abc import Generator
from typing import Type
import inspect

from pypath.internals.silver_schema import Entity as SilverEntity, Identifier, Annotation
from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    OntologyAnnotationCv,
    LicenseCV,
    UpdateCategoryCV,
    ResourceAnnotationCv,
    ResourceCv,
    CvEnum,
)
from pypath.internals import cv_terms
from ...internals.tabular_builder import (
    Annotations,
    Column,
    Entity,
    Identifiers,
)


def omnipath() -> Generator[SilverEntity]:
    """
    Yield resource metadata as an Entity record.

    Yields:
        Entity record with type CV_TERM containing OmniPath ontology metadata.
    """
    yield SilverEntity(
        type=EntityTypeCv.CV_TERM,
        identifiers=[
            Identifier(type=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=ResourceCv.OMNIPATH_ONTOLOGY),
            Identifier(type=IdentifierNamespaceCv.NAME, value='OmniPath Ontology'),
        ],
        annotations=[
            Annotation(term=ResourceAnnotationCv.LICENSE, value=str(LicenseCV.CC_BY_4_0)),
            Annotation(term=ResourceAnnotationCv.UPDATE_CATEGORY, value=str(UpdateCategoryCV.REGULAR)),
            Annotation(term=IdentifierNamespaceCv.PUBMED, value='33950214'),
            Annotation(term=ResourceAnnotationCv.URL, value='https://omnipathdb.org/'),
            Annotation(term=ResourceAnnotationCv.DESCRIPTION, value=(
                'OmniPath controlled vocabulary terms used throughout the OmniPath build system. '
                'These CV terms are organized hierarchically and include entity types, identifier '
                'namespaces, biological and experimental roles, interaction types, curation metadata, '
                'and resource metadata. Terms are based on PSI-MI standard accessions or OmniPath-specific '
                'OM accessions, with hierarchical parent relationships captured via is_a annotations.'
            )),
        ],
    )


def _extract_cv_terms() -> Generator[dict]:
    """
    Extract all CV terms from pypath.internals.cv_terms module.

    Yields:
        Dictionary with term information including accession, name, definition, and parent.
    """
    # Track parent terms that need to be added
    parent_terms = {}  # accession -> (name, definition)

    # Get all CvEnum subclasses from the cv_terms module
    for _, obj in inspect.getmembers(cv_terms):
        if inspect.isclass(obj) and issubclass(obj, CvEnum) and obj is not CvEnum:
            # Get the parent term if it exists
            parent_accession = None
            if hasattr(obj, 'parent_cv_term') and obj.parent_cv_term:
                parent_term = obj.parent_cv_term
                # parent_cv_term can be a string (MI accession) or tuple (OM accession with metadata)
                if isinstance(parent_term, tuple):
                    parent_accession = parent_term[0]
                    # If it's an OM parent term, store it for later
                    if parent_accession.startswith('OM:') and parent_accession not in parent_terms:
                        parent_name = parent_term[1] if len(parent_term) > 1 else ''
                        parent_def = parent_term[2] if len(parent_term) > 2 else ''
                        parent_terms[parent_accession] = (parent_name, parent_def)
                elif isinstance(parent_term, str):
                    parent_accession = parent_term

            # Iterate through all enum members
            for member in obj:
                # Extract member metadata
                accession = member.value

                # Only include OmniPath-specific terms (OM accessions)
                if not accession.startswith('OM:'):
                    continue

                member_name = member.name
                definition = getattr(member, 'definition', None) or ''

                # Build the term dictionary
                term_data = {
                    'accession': accession,
                    'name': member_name,
                    'definition': definition,
                    'is_a': parent_accession if parent_accession else '',
                }

                yield term_data

    # Yield parent terms that aren't already included as regular terms
    for parent_acc, (parent_name, parent_def) in parent_terms.items():
        term_data = {
            'accession': parent_acc,
            'name': parent_name,
            'definition': parent_def,
            'is_a': '',
        }
        yield term_data


def omnipath_ontology() -> Generator[SilverEntity]:
    """
    Generate OmniPath CV terms as Entity records.

    Extracts all controlled vocabulary terms from the pypath.internals.cv_terms module
    and converts them into SilverEntity records with CV_TERM type, including their
    hierarchical parent relationships.

    Yields:
        Entity records with type CV_TERM
    """
    # Define the schema mapping
    schema = Entity(
        entity_type=EntityTypeCv.CV_TERM,
        identifiers=Identifiers(
            Column('accession', cv=IdentifierNamespaceCv.CV_TERM_ACCESSION),
            Column('name', cv=IdentifierNamespaceCv.NAME),
        ),
        annotations=Annotations(
            Column('definition', cv=OntologyAnnotationCv.DEFINITION),
            Column('is_a', cv=IdentifierNamespaceCv.CV_TERM_ACCESSION),
        ),
    )

    # Extract and yield entities
    for term_data in _extract_cv_terms():
        yield schema(term_data)
