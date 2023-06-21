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

import re
import collections

import bs4

import pypath.resources.urls as urls
import pypath.share.curl as curl



def transmir_interactions():
    """
    TF-miRNA gene intractions from TransmiR.
    """

    url = urls.urls['transmir']['url']
    c = curl.Curl(
        url,
        silent = False,
        large = True,
        encoding = 'iso-8859-1',
    )

    TransmirInteraction = collections.namedtuple(
        'TransmirInteraction',
        [
            'tf_genesymbol',
            'mirna',
            'effect',
            'pubmed',
        ]
    )

    result = []

    for l in c.result:

        l = l.strip().split('\t')

        result.append(
            TransmirInteraction(
                tf_genesymbol = l[0],
                mirna = l[1],
                effect = l[4].split('(')[0],
                pubmed = l[5],
            )
        )

    return result
