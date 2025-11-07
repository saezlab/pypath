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

from pypath.internals.silver_schema import Entity, Identifier, Annotation, Membership
from pypath.internals.cv_terms import EntityTypeCv, IdentifierNamespaceCv
from ._download import (
    download_complexes,
    download_protein_families,
    download_phenotypes,
    download_stimuli,
)

__all__ = [
    'signor_complexes',
    'signor_protein_families',
    'signor_phenotypes',
    'signor_stimuli',
]


def signor_complexes() -> Generator[Entity]:
    """
    Download and parse SIGNOR complex data as Entity records.

    Yields:
        Entity records with type PROTEIN_COMPLEX, containing member proteins
    """
    for row in download_complexes():
        # Skip if we don't have at least one identifier
        if not row['complex_id']:
            continue

        # Build identifiers
        identifiers = [
            Identifier(type=IdentifierNamespaceCv.SIGNOR, value=row['complex_id'])
        ]
        if row['name']:
            identifiers.append(
                Identifier(type=IdentifierNamespaceCv.NAME, value=row['name'])
            )

        # Build member entities (proteins in the complex)
        members = []
        for component_id in row.get('components', []):
            if component_id:  # Skip empty component IDs
                member_entity = Entity(
                    source='signor',
                    type=EntityTypeCv.PROTEIN,
                    identifiers=[
                        Identifier(
                            type=IdentifierNamespaceCv.UNIPROT,
                            value=component_id
                        )
                    ],
                )
                members.append(Membership(member=member_entity))

        yield Entity(
            source='signor',
            type=EntityTypeCv.PROTEIN_COMPLEX,
            identifiers=identifiers,
            members=members if members else None,
        )


def signor_protein_families() -> Generator[Entity]:
    """
    Download and parse SIGNOR protein family data as Entity records.

    Yields:
        Entity records with type PROTEIN_FAMILY, containing member proteins
    """
    for row in download_protein_families():
        # Skip if we don't have at least one identifier
        if not row['family_id']:
            continue

        # Build identifiers
        identifiers = [
            Identifier(type=IdentifierNamespaceCv.SIGNOR, value=row['family_id'])
        ]
        if row['name']:
            identifiers.append(
                Identifier(type=IdentifierNamespaceCv.NAME, value=row['name'])
            )

        # Build member entities (proteins in the family)
        members = []
        for member_id in row.get('members', []):
            if member_id:  # Skip empty member IDs
                member_entity = Entity(
                    source='signor',
                    type=EntityTypeCv.PROTEIN,
                    identifiers=[
                        Identifier(
                            type=IdentifierNamespaceCv.UNIPROT,
                            value=member_id
                        )
                    ],
                )
                members.append(Membership(member=member_entity))

        yield Entity(
            source='signor',
            type=EntityTypeCv.PROTEIN_FAMILY,
            identifiers=identifiers,
            members=members if members else None,
        )


def signor_phenotypes() -> Generator[Entity]:
    """
    Download and parse SIGNOR phenotype data as Entity records.

    Yields:
        Entity records with type PHENOTYPE
    """
    for row in download_phenotypes():
        # Skip if we don't have at least one identifier
        if not row['phenotype_id']:
            continue

        # Build identifiers
        identifiers = [
            Identifier(type=IdentifierNamespaceCv.SIGNOR, value=row['phenotype_id'])
        ]
        if row['name']:
            identifiers.append(
                Identifier(type=IdentifierNamespaceCv.NAME, value=row['name'])
            )

        # Build annotations (description as annotation)
        annotations = []
        if row.get('description'):
            # Store description as a generic annotation
            # Note: We'd need a CV term for "description" - using NAME as synonym for now
            annotations.append(
                Annotation(
                    term=IdentifierNamespaceCv.SYNONYM,
                    value=row['description']
                )
            )

        yield Entity(
            source='signor',
            type=EntityTypeCv.PHENOTYPE,
            identifiers=identifiers,
            annotations=annotations if annotations else None,
        )


def signor_stimuli() -> Generator[Entity]:
    """
    Download and parse SIGNOR stimulus data as Entity records.

    Yields:
        Entity records with type STIMULUS
    """
    for row in download_stimuli():
        # Skip if we don't have at least one identifier
        if not row['stimulus_id']:
            continue

        # Build identifiers
        identifiers = [
            Identifier(type=IdentifierNamespaceCv.SIGNOR, value=row['stimulus_id'])
        ]
        if row['name']:
            identifiers.append(
                Identifier(type=IdentifierNamespaceCv.NAME, value=row['name'])
            )

        # Build annotations (description as annotation)
        annotations = []
        if row.get('description'):
            # Store description as a generic annotation TODO FIX: We'd need a CV term for "description" - using NAME as synonym for now
            annotations.append(
                Annotation(
                    term=IdentifierNamespaceCv.SYNONYM,
                    value=row['description']
                )
            )

        yield Entity(
            source='signor',
            type=EntityTypeCv.STIMULUS,
            identifiers=identifiers,
            annotations=annotations if annotations else None,
        )
