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

import pypath.share.curl as curl
import pypath.resources.urls as urls


AlzpathwayInteraction = collections.namedtuple(
    'AlzpathwayInteraction',
    [
        'uniprot_a',
        'uniprot_b',
        'genesymbol_a',
        'genesymbol_b',
        'pmids',
    ],
)


def alzpathway_interactions():

    url = urls.urls['alzpathway']['url']
    c = curl.Curl(url, large = True, silent = False)

    for row in c.result:

        row = row.strip().split('\t')

        yield AlzpathwayInteraction(
            uniprot_a = row[0],
            uniprot_b = row[1],
            genesymbol_a = row[4],
            genesymbol_b = row[5],
            pmids = row[8],
        )
