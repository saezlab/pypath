#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2022
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: Dénes Türei (turei.denes@gmail.com)
#           Nicolàs Palacio
#           Sebastian Lobentanzer
#           Erva Ulusoy
#           Olga Ivanova
#           Ahmet Rifaioglu
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

import collections

import pypath.resources.urls as urls
import pypath.share.curl as curl


def lncdisease_interactions():

    LncdiseaseInteraction = collections.namedtuple(
        'LncdiseaseInteraction',
        (
            'source',
            'target',
            'source_type',
            'target_type',
            'mechanism',
            'organism',
            'pmid',
        ),
    )

    url = urls.urls['lncdisease']['url_rescued']
    c = curl.Curl(url, silent = False, large = True)
    result = []

    for l in c.result:

        l = l.strip().split('\t')

        result.append(
            LncdiseaseInteraction(
                source = l[1],
                target = l[2],
                source_type = l[3].split('-')[0],
                target_type = l[3].split('-')[1] if '-' in l[3] else '',
                mechanism = l[4].lower(),
                organism = l[6].lower(),
                pmid = l[9],
            )
        )

    return result
