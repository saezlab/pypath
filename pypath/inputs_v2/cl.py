"""Cell Ontology (CL): cell types, from the basic edition without imports."""

from __future__ import annotations

from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content

from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.inputs_v2.base import (
    Dataset,
    Download,
    Resource,
    ResourceConfig,
    ontology_entity_mapper,
)
from pypath.inputs_v2.parsers.obo import iter_obo, obo_record_to_term


config = ResourceConfig(
    id=ResourceCv.CL,
    name='CL',
    full='Cell Ontology',
    url='https://obophenotype.github.io/cell-ontology/',
    # cl-basic.obo: "property_value: terms:license http://creativecommons.org/licenses/by/4.0/";
    # the OBO Foundry registry entry (obofoundry.org/ontology/cl.html) gives "CC BY 4.0",
    # same URL.
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='27377652',
    content=(Content.ONTOLOGIES,),
    summary='Cell Ontology: a structured vocabulary of cell types.',
    description=(
        'The Cell Ontology (CL) is a structured controlled vocabulary for cell types in '
        'animals. Resources cite its terms (CL:0002518) for the cell type of an experiment; '
        'this module provides their names, synonyms, definitions and hierarchy from the '
        'basic edition, which holds only CL terms and no axioms to other ontologies.'
    ),
)

terms_schema = ontology_entity_mapper(
    obo_record_to_term,
    ontology_id='cl',
    identifier_type=Namespace.CL,
    id_prefix='CL:',
)
# The basic edition's relationships are CURIEs Biolink maps (RO:0002202 develops_from,
# BFO:0000050 part_of), so no source-owned predicates are needed.


resource = Resource(
    config,
    terms=Dataset(
        download=Download(
            url='https://purl.obolibrary.org/obo/cl/cl-basic.obo',
            filename='cl-basic.obo',
            subfolder='cl',
            large=True,
        ),
        mapper=terms_schema,
        raw_parser=iter_obo,
    ),
)
