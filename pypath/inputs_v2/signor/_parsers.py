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
from ._download import (
    download_complexes,
    download_protein_families,
    download_phenotypes,
    download_stimuli,
)

SIGNOR_COMPLEX_SCHEMA = EntitySchema(
    source='signor',
    entity_type=EntityTypeCv.PROTEIN_COMPLEX,
    identifiers=Identifiers(
        Column('complex_id', cv=IdentifierNamespaceCv.SIGNOR),
        Column('name', cv=IdentifierNamespaceCv.NAME),
    ),
    members=Members(
        MemberEntities(
            entity_type=EntityTypeCv.PROTEIN,
            entity_source='signor',
            identifiers=Identifiers(
                Column('components', cv=IdentifierNamespaceCv.UNIPROT),
            ),
        )
    ),
)

SIGNOR_PROTEIN_FAMILY_SCHEMA = EntitySchema(
    source='signor',
    entity_type=EntityTypeCv.PROTEIN_FAMILY,
    identifiers=Identifiers(
        Column('family_id', cv=IdentifierNamespaceCv.SIGNOR),
        Column('name', cv=IdentifierNamespaceCv.NAME),
    ),
    members=Members(
        MemberEntities(
            entity_type=EntityTypeCv.PROTEIN,
            entity_source='signor',
            identifiers=Identifiers(
                Column('members', cv=IdentifierNamespaceCv.UNIPROT),
            ),
        )
    ),
)

SIGNOR_PHENOTYPE_SCHEMA = EntitySchema(
    source='signor',
    entity_type=EntityTypeCv.PHENOTYPE,
    identifiers=Identifiers(
        Column('phenotype_id', cv=IdentifierNamespaceCv.SIGNOR),
        Column('name', cv=IdentifierNamespaceCv.NAME),
    ),
    annotations=Annotations(
        Column('description', cv=IdentifierNamespaceCv.SYNONYM),
    ),
)

SIGNOR_STIMULUS_SCHEMA = EntitySchema(
    source='signor',
    entity_type=EntityTypeCv.STIMULUS,
    identifiers=Identifiers(
        Column('stimulus_id', cv=IdentifierNamespaceCv.SIGNOR),
        Column('name', cv=IdentifierNamespaceCv.NAME),
    ),
    annotations=Annotations(
        Column('description', cv=IdentifierNamespaceCv.SYNONYM),
    ),
)

__all__ = [
    'signor_complexes',
    'signor_protein_families',
    'signor_phenotypes',
    'signor_stimuli',
]


def signor_complexes() -> Generator[SilverEntity]:
    """
    Download and parse SIGNOR complex data as Entity records.

    Yields:
        Entity records with type PROTEIN_COMPLEX, containing member proteins
    """
    for row in download_complexes():
        entity = SIGNOR_COMPLEX_SCHEMA(row)
        if entity:
            yield entity


def signor_protein_families() -> Generator[SilverEntity]:
    """
    Download and parse SIGNOR protein family data as Entity records.

    Yields:
        Entity records with type PROTEIN_FAMILY, containing member proteins
    """
    for row in download_protein_families():
        entity = SIGNOR_PROTEIN_FAMILY_SCHEMA(row)
        if entity:
            yield entity


def signor_phenotypes() -> Generator[SilverEntity]:
    """
    Download and parse SIGNOR phenotype data as Entity records.

    Yields:
        Entity records with type PHENOTYPE
    """
    for row in download_phenotypes():
        entity = SIGNOR_PHENOTYPE_SCHEMA(row)
        if entity:
            yield entity


def signor_stimuli() -> Generator[SilverEntity]:
    """
    Download and parse SIGNOR stimulus data as Entity records.

    Yields:
        Entity records with type STIMULUS
    """
    for row in download_stimuli():
        entity = SIGNOR_STIMULUS_SCHEMA(row)
        if entity:
            yield entity
