"""
Parse SLC tables data and emit Entity records.

This module converts SLC tables reference data into Entity records using the
declarative schema pattern.
"""

import requests
from functools import reduce

from bs4 import BeautifulSoup

from pypath.inputs_v2.parsers.base import iter_csv
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
    IdentifierNamespaceCv,
    LicenseCV,
    UpdateCategoryCV,
    ResourceCv,
    MoleculeAnnotationsCv,
)

# XXX: Do we need just base tables or also the info inside each href of each SLC

# =================================== SET-UP ===================================

BASE_URL = 'https://slc.bioparadigms.org/'


def chunk_this(L, n):
    '''
    For a given list *L*, returns another list of *n*-sized chunks from
    it (in the same order).

    * Arguments:
        - *L* [list]: The list to be sliced into sublists of the
          definded size.
        - *n* [int]: The size of the chunks.

    * Returns:
        - [list]: List of *n*-sized chunks from *L*. **NOTE:** If the
          number of items in *L* is not divisible by *n*, the last
          element returned will have an inferior size.

    * Examples:
        >>> L = range(6)
        >>> chunk_this(L, 2)
        [[0, 1], [2, 3], [4, 5]]
        >>> chunk_this(L, 4)
        [[0, 1, 2, 3], [4, 5]]
    '''
    L = list(L)

    return [L[i:i + n] for i in range(0, len(L), n)]


def parser(opener, **_kwargs):

    txt = reduce(str.__add__, list(opener.result))

    soup = BeautifulSoup(txt, 'html.parser')

    tables = soup.find_all('table')
    slcfams = [t.text for t in soup.find_all('span', attrs={'class': 'slcname'})]
    raw_tables = soup.find_all('table')

    headers = ['SLC_family'] + [
        t.text for t in raw_tables[0].find_all('td', attrs={'class': 'tbl_head'})
    ]

    rows = []

    for i, slcf in enumerate(slcfams):

        tbl = raw_tables[i]
        rows += [[slcf] + chunk for chunk in chunk_this(
            [t.text for t in tbl.find_all('td', attrs={'class': 'tbl_cell'})],
            len(headers) - 1
        )]

    yield from iter(headers + rows)


config = ResourceConfig(
    id=ResourceCv.SLC_TABLES,
    name='SLC Tables',
    url='https://slc.bioparadigms.org/',
    license=LicenseCV.UNSPECIFIED,
    update_category=UpdateCategoryCV.STATIC, # I guess?
    pubmed=''# Dunno, gotta search coz not mentioned on website X(
)

# ================================== DOWNLOAD ==================================

download = Download(
    url=BASE_URL,
    filename='slctables.html',
    subfolder='slctables',
    large=True,
    ext='.html',
    default_mode='r',
)


# ================================= REFERENCE ==================================

# SLC_family                        SLC1
# SLC name                          SLC1A1
# Protein name                      EAAC1, EAAT3
# Aliases                           System X-AG
# Transport type*                   C / Na+, H+, K+
# Substrates                        L-Glu, D/L-Asp
# Tissue and cellular expression    brain (neurons), intestine, kidney, ...