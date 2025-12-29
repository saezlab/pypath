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

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    OntologyAnnotationCv,
    LicenseCV,
    UpdateCategoryCV,
    ResourceCv,
)
from ...internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
)
from .shared import iter_obo_terms


# PSI-MI OBO URL
PSI_MI_URL = "https://raw.githubusercontent.com/HUPO-PSI/psi-mi-CV/master/psi-mi.obo"

config = ResourceConfig(
    id=ResourceCv.PSI_MI,
    name='PSI-MI Controlled Vocabulary',
    url='https://github.com/HUPO-PSI/psi-mi-CV',
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='17925023',
    description=(
        'The PSI-MI (Proteomics Standards Initiative - Molecular Interactions) '
        'controlled vocabulary provides standardized terms for describing '
        'molecular interactions. It includes terms for interaction types, '
        'detection methods, participant roles, and experimental features.'
    ),
)

f = FieldConfig()
schema = EntityBuilder(
    entity_type=EntityTypeCv.CV_TERM,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=f('accession')),
        CV(term=IdentifierNamespaceCv.NAME, value=f('name')),
        CV(term=IdentifierNamespaceCv.SYNONYM, value=f('synonyms', delimiter=';')),
        CV(term=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=f('alt_ids', delimiter=';')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=OntologyAnnotationCv.DEFINITION, value=f('definition')),
        CV(term=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=f('is_a', delimiter=';')),
        CV(term=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=f('xrefs', delimiter=';')),
        CV(term=OntologyAnnotationCv.COMMENT, value=f('comment')),
        CV(term=OntologyAnnotationCv.IS_OBSOLETE, value=f('is_obsolete')),
    ),
)

resource = Resource(
    config,
    psi_mi_ontology=Dataset(
        download=Download(
            url=PSI_MI_URL,
            filename='psi-mi.obo',
            subfolder='psi_mi',
        ),
        mapper=schema,
        raw_parser=iter_obo_terms,
    ),
)

psi_mi_ontology = resource.psi_mi_ontology
