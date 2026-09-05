"""Gene Ontology term resource."""

from __future__ import annotations

from biolink_model.datamodel.model import slots
from omnipath_core.naming import Namespace

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
    id=ResourceCv.GENE_ONTOLOGY,
    name='Gene Ontology',
    url='https://geneontology.org/',
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='41413728',
    primary_category='ontologies',
    resource_kind='ontology',
    description='The Gene Ontology provides controlled vocabulary terms for gene product function, process, and cellular location.',
)

_GENE_ONTOLOGY_ID = 'gene_ontology'
terms_schema = ontology_entity_mapper(
    obo_record_to_term,
    ontology_id=_GENE_ONTOLOGY_ID,
    identifier_type=Namespace.GO,
    # GO explicitly uses these BFO/symbolic predicates for mereology.
    relationship_predicates={
        'part_of': slots.part_of, 'BFO:0000050': slots.part_of,
        'has_part': slots.has_part, 'BFO:0000051': slots.has_part,
    },
)


resource = Resource(
    config,
    terms=Dataset(
        download=Download(
            url='https://purl.obolibrary.org/obo/go.obo',
            filename='go.obo',
            subfolder='go',
            large=True,
        ),
        mapper=terms_schema,
        raw_parser=iter_obo,
    ),
)
