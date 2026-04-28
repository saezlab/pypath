"""Gene Ontology OBO artifact resource."""

from __future__ import annotations

from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.inputs_v2.base import (
    ArtifactDataset,
    Download,
    Resource,
    ResourceConfig,
    read_opener_text,
)


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
    ontology=ArtifactDataset(
        download=Download(
            url='https://purl.obolibrary.org/obo/go.obo',
            filename='go.obo',
            subfolder='go',
            large=True,
        ),
        renderer=read_opener_text,
        extension='obo',
        file_stem='go',
        kind='ontology',
    ),
)
