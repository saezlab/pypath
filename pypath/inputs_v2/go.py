"""Gene Ontology OBO artifact resource."""

from __future__ import annotations

from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.inputs_v2.base import Download, OntologyDataset, Resource, ResourceConfig
from pypath.inputs_v2.parsers.obo import iter_obo, obo_record_to_term


config = ResourceConfig(
    id=ResourceCv.GENE_ONTOLOGY,
    name='Gene Ontology',
    url='https://geneontology.org/',
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    primary_category='ontologies',
    resource_kind='ontology',
    description='The Gene Ontology provides controlled vocabulary terms for gene product function, process, and cellular location.',
)


resource = Resource(
    config,
    ontology=OntologyDataset(
        download=Download(
            url='https://purl.obolibrary.org/obo/go.obo',
            filename='go.obo',
            subfolder='go',
            large=True,
        ),
        mapper=obo_record_to_term,
        raw_parser=iter_obo,
        ontology_id='gene_ontology',
        extension='obo',
        file_stem='go',
    ),
)
