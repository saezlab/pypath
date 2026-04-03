"""
Parse RaMP-DB data and emit Entity records.

This module converts annotations of lipids and metabolites into Entity records
using the declarative schema pattern.
"""

import os
import re
import requests

from bs4 import BeautifulSoup

from pypath.inputs_v2.base import ResourceConfig
from pypath.internals.cv_terms import ResourceCv, LicenseCV, UpdateCategoryCV
from pypath.inputs_v2.base import Download


def get_ramp_latest_ver(branch='main'):

    url = f'https://github.com/ncats/RaMP-DB/raw/refs/heads/{branch}/db/'
    res = requests.get(url)
    soup = BeautifulSoup(res.text, 'html.parser')
    files = sorted({
        f.text for f in soup.find_all(title=re.compile("\\.sqlite.gz$"))
    })

    return url + files[-1]


config = ResourceConfig(
    id=ResourceCv.RAMPDB,
    name='RaMP-DB',
    url='https://rampdb.nh.gov/',
    license=LicenseCV.GPL_2_0,
    update_category=UpdateCategoryCV.IRREGULAR,
    pubmed='36373969',
    primary_category='pathways',
    description=(
        'RaMP-DB (Relational database of Metabolomic Pathways) is a '
        'multi-sourced integrated database with comprehensive annotations on '
        'biological pathways, structure/chemistry, disease and ontology '
        'annotations for genes, proteins, and metabolites.'
    )
)

url = get_ramp_latest_ver()

download = Download(
    url=url,
    filename=os.path.basename(url),
    subfolder='ramp',
    large=True,
    ext=url.split('.')[-1],
)
