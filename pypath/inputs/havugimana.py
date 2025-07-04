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

import pypath.inputs.cell as cell_input
import pypath.inputs.common as inputs_common
import pypath.resources.urls as urls
import pypath.internals.intera as intera


def get_havugimana():
    """
    Downloads data from
    Supplement Table S3/1 from Havugimana 2012
    Cell. 150(5): 1068–1081.
    """

    # Use url_rescued (local file) if available, otherwise try online download
    url = urls.urls['havugimana'].get('url_rescued', urls.urls['havugimana']['url'])
    if 'url_rescued' in urls.urls['havugimana']:
        path = url  # Use local file path directly
        table = inputs_common.read_xls(path)
    else:
        # Fallback to original method if rescued URL not available
        article_url = urls.urls['havugimana']['article']
        path = cell_input.cell_supplementary(url, article_url)
        table = inputs_common.read_xls(path)
        
    return table[3:]  # Skip header rows


def havugimana_complexes():
    """
    Retrieves complexes from
    Supplement Table S3/1 from Havugimana 2012
    Cell. 150(5): 1068–1081.
    """

    complexes = {}

    for rec in get_havugimana():

        cplex = intera.Complex(
            components = rec[2].split(','),
            sources = 'Havugimana2012',
            ids = rec[0],
        )

        complexes[cplex.__str__()] = cplex

    return complexes
