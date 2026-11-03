"""Human Disease Ontology (DOID): human diseases, classified by etiology."""

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
    id=ResourceCv.DOID,
    name='DOID',
    full='Human Disease Ontology',
    url='https://disease-ontology.org/',
    # doid.obo: "property_value: terms:license
    # https://creativecommons.org/publicdomain/zero/1.0/" and "remark: The Disease Ontology
    # content is available via the Creative Commons Public Domain Dedication CC0 1.0 Universal
    # license"; the OBO Foundry registry entry (obofoundry.org/ontology/doid.html) gives
    # "CC0 1.0", https://creativecommons.org/publicdomain/zero/1.0/.
    license=LicenseCV.CC0_1_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='34755882',
    content=(Content.ONTOLOGIES,),
    summary='Human Disease Ontology: human diseases classified by etiology, with synonyms and cross-references to disease vocabularies.',
    description=(
        'The Human Disease Ontology (DOID), maintained by the Institute for Genome Sciences '
        'of the University of Maryland School of Medicine, is a standardized ontology of '
        'human disease: infectious, genetic, cancer, metabolic and other diseases, '
        'classified by etiology, with synonyms and cross-references to MeSH, ICD, OMIM, '
        'NCI Thesaurus, SNOMED CT and UMLS. Resources cite its terms (DOID:3702) for the '
        'disease a cell line derives from; this module provides their names, synonyms, '
        'definitions and hierarchy.'
    ),
)

# doid.obo holds only DOID terms, the classified hierarchy (asserted parents and those a
# reasoner infers from the logical definitions); the prefix check keeps it so were a
# foreign term ever merged in.
terms_schema = ontology_entity_mapper(
    obo_record_to_term,
    ontology_id='doid',
    identifier_type=Namespace.DOID,
    id_prefix='DOID:',
)


resource = Resource(
    config,
    terms=Dataset(
        download=Download(
            url='https://purl.obolibrary.org/obo/doid.obo',
            filename='doid.obo',
            subfolder='doid',
            large=True,
        ),
        mapper=terms_schema,
        raw_parser=iter_obo,
    ),
)
