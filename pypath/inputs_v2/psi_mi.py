"""PSI-MI controlled vocabulary OBO artifact resource."""

from __future__ import annotations

import re

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
    id=ResourceCv.PSI_MI,
    name='PSI-MI Controlled Vocabulary',
    url='https://github.com/HUPO-PSI/psi-mi-CV',
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='14755292',
    primary_category='ontologies',
    resource_kind='ontology',
    description='The PSI-MI controlled vocabulary provides molecular interaction terms and metadata.',
)


def _fix_obo_text(content: str) -> str:
    content = re.sub(r'^date: \d{2}:\d{2}:\d{4}.*\n', '', content, flags=re.MULTILINE)
    return re.sub(r'^creation_date:.*\n', '', content, flags=re.MULTILINE)


def _iter_fixed_obo(opener, **kwargs):
    yield from iter_obo(opener, preprocess_text=_fix_obo_text, **kwargs)


terms_schema = ontology_entity_mapper(obo_record_to_term, ontology_id='psi-mi')


resource = Resource(
    config,
    terms=Dataset(
        download=Download(
            url='https://raw.githubusercontent.com/HUPO-PSI/psi-mi-CV/master/psi-mi.obo',
            filename='psi-mi.obo',
            subfolder='psi_mi',
            large=True,
        ),
        mapper=terms_schema,
        raw_parser=_iter_fixed_obo,
    ),
)
