"""PSI-MI controlled vocabulary OBO artifact resource."""

from __future__ import annotations

import re

from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.inputs_v2.base import (
    ArtifactDataset,
    Download,
    Resource,
    ResourceConfig,
    read_opener_text,
)


config = ResourceConfig(
    id=ResourceCv.PSI_MI,
    name='PSI-MI Controlled Vocabulary',
    url='https://github.com/HUPO-PSI/psi-mi-CV',
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    primary_category='ontologies',
    resource_kind='ontology',
    description='The PSI-MI controlled vocabulary provides molecular interaction terms and metadata.',
)


def _read_and_fix_obo(opener, **_kwargs) -> str:
    content = read_opener_text(opener)
    content = re.sub(r'^date: \d{2}:\d{2}:\d{4}.*\n', '', content, flags=re.MULTILINE)
    return re.sub(r'^creation_date:.*\n', '', content, flags=re.MULTILINE)


resource = Resource(
    config,
    ontology=ArtifactDataset(
        download=Download(
            url='https://raw.githubusercontent.com/HUPO-PSI/psi-mi-CV/master/psi-mi.obo',
            filename='psi-mi.obo',
            subfolder='psi_mi',
            large=True,
        ),
        renderer=_read_and_fix_obo,
        extension='obo',
        file_stem='psi_mi',
        kind='ontology',
    ),
)
