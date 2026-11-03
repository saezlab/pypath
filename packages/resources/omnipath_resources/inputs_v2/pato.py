"""Phenotype And Trait Ontology (PATO): the qualities that phenotypes are described with."""

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
    id=ResourceCv.PATO,
    name='PATO',
    full='Phenotype And Trait Ontology',
    url='https://github.com/pato-ontology/pato/',
    # pato.obo: "property_value: terms:license https://creativecommons.org/licenses/by/3.0/";
    # the OBO Foundry registry entry (obofoundry.org/ontology/pato.html) gives "CC BY 3.0",
    # http://creativecommons.org/licenses/by/3.0/.
    license=LicenseCV.CC_BY_3_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='28387809',
    content=(Content.ONTOLOGIES,),
    summary='Phenotype And Trait Ontology: phenotypic qualities (size, shape, color, viability, increased, decreased) used to describe phenotypes.',
    description=(
        'The Phenotype And Trait Ontology (PATO) is an ontology of phenotypic qualities: the '
        'properties, attributes and characteristics of an organism or its parts, such as '
        'viability, size or color, and their values (increased, decreased, absent). Other '
        'ontologies compose them with anatomy or process terms to describe a phenotype. '
        'Resources cite its terms (PATO:0000169) for the phenotype of a genetic interaction '
        'or the environment of an experiment; this module provides their names, synonyms, '
        'definitions and hierarchy.'
    ),
)

terms_schema = ontology_entity_mapper(
    obo_record_to_term,
    ontology_id='pato',
    identifier_type=Namespace.PATO,
    id_prefix='PATO:',
    # PATO writes parthood and oppositeness as symbolic relationship types Biolink does not see.
    relationship_predicates={
        'has_part': slots.has_part,
        'is_opposite_of': slots.opposite_of,
    },
)


resource = Resource(
    config,
    terms=Dataset(
        download=Download(
            # pato.obo holds the PATO terms alone; pato-base.obo adds a stub of a Cell Ontology term.
            url='https://purl.obolibrary.org/obo/pato.obo',
            filename='pato.obo',
            subfolder='pato',
            large=True,
        ),
        mapper=terms_schema,
        raw_parser=iter_obo,
    ),
)
