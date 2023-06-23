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

import itertools
import collections

import pypath.resources.urls as urls
import pypath.share.curl as curl


def mir2disease_interactions():
    """
    Retrieves literature curated miRNA-target gene interactions from the
    miR2Disease database (http://www.mir2disease.org/).
    """

    Mir2diseaseInteraction = collections.namedtuple(
        'Mir2diseaseInteraction',
        (
            'mirna',
            'target_genesymbol',
            'year',
            'sentence',
        ),
    )

    url = urls.urls['mir2dis']['url_rescued']
    c = curl.Curl(url, silent = True, large = True, encoding = 'iso-8859-1')

    return [
        Mir2diseaseInteraction(
            *l.strip('\r\n\t "').split('\t')
        )
        for l in itertools.islice(c.result, 3, None)
        if l
    ]
