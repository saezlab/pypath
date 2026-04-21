"""
Parse MACdb data and emit Entity records.

This module converts annotations of metabolites-cancer associations into Entity
records using the declarative schema pattern.
"""

import os
import re
import requests
from functools import partial
from pathlib import Path

from bs4 import BeautifulSoup

from pypath.inputs_v2.parsers.base import iter_tsv
from pypath.inputs_v2.base import ResourceConfig, Download, Resource, Dataset
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
)
from pypath.internals.cv_terms import (
    EntityTypeCv,
    BiologicalRoleCv,
    IdentifierNamespaceCv,
    LicenseCV,
    UpdateCategoryCV,
    ResourceCv,
    MoleculeSubtypeCv,
    MoleculeAnnotationsCv,
    AssayAnnotationsCv,
)

# =================================== SET-UP ===================================

BASE_URL = 'https://ngdc.cncb.ac.cn/macdb/api/get_download_file?file_type=%s'
TABLES = ['metabolite', 'trait', 'study', 'publication']

config = ResourceConfig(
    id=ResourceCv.MACDB,
    name='MACdb',
    url='https://ngdc.cncb.ac.cn/macdb/',
    license=LicenseCV.ACADEMIC_FREE,
    update_category=UpdateCategoryCV.IRREGULAR, # Static?
    pubmed='37027007',
    primary_category='metabolism',
    description=(
        'MACdb is a comprehensive knowledgebase based totally on manual '
        'curation of cancer metabolic literature, providing high quality '
        'metabolite associations, online visualizing networks, and other '
        'effective tools for cancer or metabonomics researches.'
    ),
)

# ================================== DOWNLOAD ==================================

download = dict()

for t in TABLES:

    download[t] = Download(
            url=BASE_URL % t,
            filename='downloads.%s.txt',
            subfolder='macdb',
            large=True,
            ext='.txt',
            default_mode='r',
        )

# =================================== SCHEMA ===================================

f = FieldConfig(
    extract={},
    map={},
    transform={},
)

metabolite_schema = EntityBuilder(

)

trait_schema = EntityBuilder()

study_schema = EntityBuilder()

publication_schema = EntityBuilder()

# ================================= RESOURCE ===================================
