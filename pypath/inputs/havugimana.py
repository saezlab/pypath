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

import collections
import itertools

import pypath.share.curl as curl
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

    supp_url = urls.urls['havugimana']['url']
    article_url = urls.urls['havugimana']['article']
    path = cell_input.cell_supplementary(supp_url, article_url)
    table = inputs_common.read_xls(path)

    return table[3:]


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
