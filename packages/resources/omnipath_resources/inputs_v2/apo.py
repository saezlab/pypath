"""Ascomycete Phenotype Ontology (APO): the yeast phenotypes of SGD and their qualifiers."""

from __future__ import annotations

from typing import Any

from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content

from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.ontology_schema import OntologyTerm
from pypath.inputs_v2.base import (
    Dataset,
    Download,
    Resource,
    ResourceConfig,
    ontology_entity_mapper,
)
from pypath.inputs_v2.parsers.obo import iter_obo, obo_record_to_term


config = ResourceConfig(
    id=ResourceCv.APO,
    name='APO',
    full='Ascomycete Phenotype Ontology',
    url='https://www.yeastgenome.org/',
    # apo.obo: "property_value: http://purl.org/dc/terms/license
    # http://creativecommons.org/licenses/by/4.0/"; the OBO Foundry registry entry
    # (obofoundry.org/ontology/apo.html) gives "CC BY 4.0", https://creativecommons.org/licenses/by/4.0/.
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='20157474',
    content=(Content.ONTOLOGIES,),
    summary='Ascomycete Phenotype Ontology: the phenotypes observed in yeast mutants, with their qualifiers, mutant types and experiment types.',
    description=(
        'The Ascomycete Phenotype Ontology (APO), maintained by the Saccharomyces Genome '
        'Database, is a structured controlled vocabulary for the phenotypes of Ascomycete '
        'fungi. It holds observables (colony size, vegetative growth), qualifiers (abnormal, '
        'increased), mutant types and experiment types. Resources cite its terms (APO:0000063) '
        'for the phenotype of a genetic interaction; this module provides their names, '
        'synonyms, definitions and hierarchy.'
    ),
)


def apo_record_to_term(row: dict[str, Any]) -> OntologyTerm | None:
    """The term of a record, without the alt_ids that are not accessions.

    Two terms of the file carry ``alt_id: curators`` (a curator's name, not an id); it would
    be one identifier shared by both of them. The ``YPO:`` alt_ids (the former yeast
    phenotype ontology's ids of the same terms) are kept.
    """
    alt_ids = [value for value in row.get('alt_ids') or [] if ':' in value]
    return obo_record_to_term({**row, 'alt_ids': alt_ids})


terms_schema = ontology_entity_mapper(
    apo_record_to_term,
    ontology_id='apo',
    identifier_type=Namespace.APO,
    id_prefix='APO:',
)


resource = Resource(
    config,
    terms=Dataset(
        download=Download(
            url='https://purl.obolibrary.org/obo/apo.obo',
            filename='apo.obo',
            subfolder='apo',
            large=True,
        ),
        mapper=terms_schema,
        raw_parser=iter_obo,
    ),
)
