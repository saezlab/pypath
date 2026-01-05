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


# Gene Ontology OBO URL
GENE_ONTOLOGY_URL = "https://current.geneontology.org/ontology/go.obo"

config = ResourceConfig(
    id=ResourceCv.GENE_ONTOLOGY,
    name='Gene Ontology',
    url='http://geneontology.org/',
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='33290552',
    description=(
        'The Gene Ontology (GO) provides a comprehensive framework for '
        'describing gene functions across all organisms. GO covers three domains: '
        'biological process, molecular function, and cellular component. '
        'GO terms are organized in a hierarchical structure with rich relationships.'
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
        CV(term=OntologyAnnotationCv.IS_OBSOLETE, value=f('is_obsolete')),
    ),
)

resource = Resource(
    config,
    gene_ontology=Dataset(
        download=Download(
            url=GENE_ONTOLOGY_URL,
            filename='go.obo',
            subfolder='gene_ontology',
        ),
        mapper=schema,
        raw_parser=iter_obo_terms,
    ),
)
