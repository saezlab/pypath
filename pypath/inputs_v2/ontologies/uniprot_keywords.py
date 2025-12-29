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

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    OntologyAnnotationCv,
    ResourceCv,
    LicenseCV,
    UpdateCategoryCV,
)
from ...internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
)
from .shared import iter_obo_terms


# UniProt Keywords OBO URL
UNIPROT_KEYWORDS_URL = "https://rest.uniprot.org/keywords/stream?format=obo&query=%28*%29"

config = ResourceConfig(
    id=ResourceCv.UNIPROT_KEYWORDS,
    name='UniProt Keywords',
    url='https://www.uniprot.org/keywords/',
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='33237286',
    description=(
        'UniProt Keywords provide a controlled vocabulary for '
        'summarizing protein properties, including biological processes, '
        'molecular functions, cellular components, protein families, and more. '
        'Keywords facilitate consistent annotation and efficient searching.'
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
        CV(term=OntologyAnnotationCv.COMMENT, value=f('comment')),
        CV(term=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=f('is_a', delimiter=';')),
        CV(term=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=f('xrefs', delimiter=';')),
    ),
)

resource = Resource(
    config,
    uniprot_keywords_ontology=Dataset(
        download=Download(
            url=UNIPROT_KEYWORDS_URL,
            filename='uniprot_keywords.obo',
            subfolder='uniprot_keywords',
        ),
        mapper=schema,
        raw_parser=iter_obo_terms,
    ),
)

uniprot_keywords_ontology = resource.uniprot_keywords_ontology
