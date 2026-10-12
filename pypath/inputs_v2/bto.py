"""BRENDA Tissue Ontology (BTO): tissues, cell lines, cell types and cell cultures."""

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
    id=ResourceCv.BTO,
    name='BTO',
    full='BRENDA Tissue Ontology',
    url='https://www.brenda-enzymes.org/ontology.php?ontology_id=3',
    # bto.obo: "property_value: http://purl.org/dc/terms/license
    # https://creativecommons.org/licenses/by/4.0/"; the OBO Foundry registry entry
    # (obofoundry.org/ontology/bto.html) gives "CC BY 4.0", same URL.
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='21030441',
    content=(Content.ONTOLOGIES,),
    summary='BRENDA Tissue Ontology: tissues, cell lines, cell types and cell cultures that are enzyme sources.',
    description=(
        'The BRENDA Tissue Ontology (BTO) is a structured controlled vocabulary for the '
        'source of an enzyme, comprising tissues, cell lines, cell types and cell cultures '
        'of uni- and multicellular organisms. Resources cite its terms (BTO:0000007) for the '
        'cell line or tissue of an experiment; this module provides their names, synonyms, '
        'definitions and hierarchy.'
    ),
)

terms_schema = ontology_entity_mapper(
    obo_record_to_term,
    ontology_id='bto',
    identifier_type=Namespace.BTO,
    id_prefix='BTO:',
    # BTO writes part of and development as symbolic relationship types Biolink does not see.
    relationship_predicates={
        'part_of': slots.part_of,
        'develops_from': slots.develops_from,
    },
)


resource = Resource(
    config,
    terms=Dataset(
        download=Download(
            url='https://purl.obolibrary.org/obo/bto.obo',
            filename='bto.obo',
            subfolder='bto',
            large=True,
        ),
        mapper=terms_schema,
        raw_parser=iter_obo,
    ),
)
