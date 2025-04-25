#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright 2014-2023
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: see the file `README.rst`
#  Contact: Dénes Türei (turei.denes@gmail.com)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      https://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: https://pypath.omnipathdb.org/
#

from typing import Literal
from collections.abc import Generator

import json

from pypath.share import curl
from pypath.resources import urls

DATA = Literal[
    "target",
    "assay",
    "molecule",
    "activity",
    "document",
    "drug_indication",
    "mechanism",
]

def json_pages(data_type: DATA,
                   max_pages: int | None = None) -> Generator[dict]:
    """
    Retrieves data from ChEMBL.

    This function will retrieve the specified data from ChEMBL
    and yield the data as a JSON object.

    Args:
        data_type (DATA): The type of data to retrieve.
        max_pages (int): The maximum number of pages to retrieve.

    Yields:
        dict: The JSON object of the retrieved data.
    """

    page_dict: dict = {}
    page_count = 0

    while True:

        if not page_dict:

            url_base = urls.urls['chembl']['url']
            url_path = urls.urls['chembl'][data_type]
            url = f"{url_base}{url_path}"

        elif page_dict['page_meta']['next']:

            url_base = urls.urls['chembl']['url']
            url_path = page_dict['page_meta']['next']
            url = f"{url_base}{url_path}"

        else:
            break

        c = curl.Curl(url, large=True, silent=False)
        with open(c.fileobj.name, mode="r", encoding='utf-8') as read_file:
            page_dict = json.load(read_file) # load the json object

            # remove unwanted page meta data
            for key in page_dict.keys():
                if key == 'page_meta':
                    continue
                else:
                    data_dict = page_dict[key]

                    # extract data from each page
                    for data in data_dict:
                        yield data

        # allows user to specify maximum number of pages
        page_count += 1
        if max_pages and page_count >= max_pages:
            break
