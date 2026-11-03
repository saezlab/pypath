"""C. elegans Phenotype Ontology (WBPhenotype): the phenotypes of worms."""

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
    id=ResourceCv.WBPHENOTYPE,
    name='WBPhenotype',
    full='C. elegans Phenotype Ontology',
    url='https://github.com/obophenotype/c-elegans-phenotype-ontology',
    # wbphenotype.obo: "property_value: terms:license http://creativecommons.org/licenses/by/4.0/";
    # the OBO Foundry registry entry (obofoundry.org/ontology/wbphenotype.html) gives "CC BY 4.0",
    # http://creativecommons.org/licenses/by/4.0/.
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='21261995',
    content=(Content.ONTOLOGIES,),
    summary='C. elegans Phenotype Ontology: the observable phenotypes of Caenorhabditis elegans and other nematodes.',
    description=(
        'The C. elegans Phenotype Ontology (WBPhenotype), maintained by WormBase, is a '
        'structured controlled vocabulary of the phenotypes of Caenorhabditis elegans and '
        'other nematodes: variants of anatomy, physiology, behavior and development. '
        'Resources cite its terms (WBPhenotype:0000050) for the phenotype of a genetic '
        'interaction; this module provides their names, synonyms, definitions and hierarchy. '
        'The prefix is written in mixed case, as the ontology and its users write it.'
    ),
)

terms_schema = ontology_entity_mapper(
    obo_record_to_term,
    ontology_id='wbphenotype',
    identifier_type=Namespace.WBPHENOTYPE,
    id_prefix='WBPhenotype:',
)


resource = Resource(
    config,
    terms=Dataset(
        download=Download(
            url='https://purl.obolibrary.org/obo/wbphenotype.obo',
            filename='wbphenotype.obo',
            subfolder='wbphenotype',
            large=True,
        ),
        mapper=terms_schema,
        raw_parser=iter_obo,
    ),
)
