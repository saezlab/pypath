"""Mammalian Phenotype Ontology (MP): the phenotypes of mice and other mammals."""

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
    id=ResourceCv.MP,
    name='MP',
    full='Mammalian Phenotype Ontology',
    url='https://www.informatics.jax.org/vocab/mp_ontology',
    # mp.obo: "property_value: terms:license https://creativecommons.org/licenses/by/4.0/";
    # the OBO Foundry registry entry (obofoundry.org/ontology/mp.html) gives "CC BY 4.0",
    # http://creativecommons.org/licenses/by/4.0/.
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='40963406',
    content=(Content.ONTOLOGIES,),
    summary='Mammalian Phenotype Ontology: the standard terms for the phenotypes of mice and other mammals.',
    description=(
        'The Mammalian Phenotype Ontology (MP), maintained by the Mouse Genome Informatics '
        'group at the Jackson Laboratory, is a structured vocabulary of the phenotypes '
        'observed in mammals, used to annotate mouse mutants and the results of '
        'high-throughput phenotyping. Resources cite its terms (MP:0001262) for the '
        'phenotype of a genetic interaction; this module provides their names, synonyms, '
        'definitions and hierarchy.'
    ),
)

# mp.obo, the release MGI publishes at the PURL, holds only MP terms and writes no
# relationship but is_a (mp-base.obo adds five 'has modifier' relationships to PATO terms).
terms_schema = ontology_entity_mapper(
    obo_record_to_term,
    ontology_id='mp',
    identifier_type=Namespace.MP,
    id_prefix='MP:',
)


resource = Resource(
    config,
    terms=Dataset(
        download=Download(
            url='https://purl.obolibrary.org/obo/mp.obo',
            filename='mp.obo',
            subfolder='mp',
            large=True,
        ),
        mapper=terms_schema,
        raw_parser=iter_obo,
    ),
)
