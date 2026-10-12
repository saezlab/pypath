"""Uberon multi-species anatomy ontology, from the basic edition without imports."""

from __future__ import annotations

from biolink_model.datamodel.model import slots
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
    id=ResourceCv.UBERON,
    name='Uberon',
    full='Uberon multi-species anatomy ontology',
    url='https://obophenotype.github.io/uberon/',
    # uberon-basic.obo: "property_value: dcterms-license
    # http://creativecommons.org/licenses/by/3.0/"; the OBO Foundry registry entry
    # (obofoundry.org/ontology/uberon.html) gives "CC BY 3.0", same URL.
    license=LicenseCV.CC_BY_3_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='22293552',
    content=(Content.ONTOLOGIES,),
    summary='Uberon: an integrated cross-species anatomy ontology of animal structures.',
    description=(
        'Uberon is an integrated cross-species anatomy ontology covering animals and '
        'bridging species-specific anatomy ontologies. Resources cite its terms '
        '(UBERON:0002113) for the tissue or organ of an experiment; this module provides '
        'their names, synonyms, definitions and hierarchy from the basic edition, which '
        'holds only Uberon terms and no imported ones.'
    ),
)

terms_schema = ontology_entity_mapper(
    obo_record_to_term,
    ontology_id='uberon',
    identifier_type=Namespace.UBERON,
    id_prefix='UBERON:',
    # The basic edition writes part of and development as symbolic relationship types.
    relationship_predicates={
        'part_of': slots.part_of,
        'develops_from': slots.develops_from,
    },
)


resource = Resource(
    config,
    terms=Dataset(
        download=Download(
            url='https://purl.obolibrary.org/obo/uberon/uberon-basic.obo',
            filename='uberon-basic.obo',
            subfolder='uberon',
            large=True,
        ),
        mapper=terms_schema,
        raw_parser=iter_obo,
    ),
)
