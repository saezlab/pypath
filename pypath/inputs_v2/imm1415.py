"""
Parse iMM1415 data and emit Entity records.

This module converts information from the BiGG Model iMM1415 into Entity
records using the declarative schema pattern.
"""

import re
from collections import defaultdict

from pypath.inputs_v2.base import (
    ResourceConfig,
    Download,
    Resource,
    Dataset,
)
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
)
from pypath.internals.cv_terms import (
    EntityTypeCv,
    LigandTypeCv,
    InteractionParameterCv,
    IdentifierNamespaceCv,
    LicenseCV,
    UpdateCategoryCV,
    ResourceCv,
)

# =================================== SET-UP ===================================

URL = 'http://bigg.ucsd.edu/static/models/iMM1415.json'

config = ResourceConfig(
    id=ResourceCv.IMM1415,
    name='iMM1415',
    url='http://bigg.ucsd.edu/models/iMM1415',
    license=LicenseCV.BIGG,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='20959003',
    primary_category='metabolism',
    description=(
        'A high-quality mouse genome-scale metabolic reconstruction, iMM1415'
        '(Mus Musculus, 1415 genes)'
    ),
)

# ================================== DOWNLOAD ==================================

download = Download(
    url=URL,
    filename='iMM1415.json',
    subfolder='bigg',
    large=True,
    ext='json',
    default_mode='r',
)

# =================================== SCHEMA ===================================

f = FieldConfig(
    extract={},
    map={},
    transform={},
)

#schema = EntityBuilder()

# ================================= RESOURCE ===================================

#resource = Resource(
#    config=config,
#    data=Dataset(
#        download=download,
#        mapper=schema,
#        raw_parser=parser,
#    ),
#)

# ================================= REFERENCE ==================================
