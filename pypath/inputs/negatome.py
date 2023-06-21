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

import pypath.resources.urls as urls
import pypath.share.curl as curl


def negatome_interactions():
    """
    Literature curated non-interacting protein pairs from the Negatome
    database. These interactions do not exist to best of our knowledge,
    the literature references point to papers with experiments testing
    for the interaction but finding negative outcome.
    """

    NegatomeInteraction = collections.namedtuple(
        'NegatomeInteraction',
        (
            'uniprot_a',
            'uniprot_b',
            'pmid',
            'method',
        ),
    )

    url = urls.urls['negatome']['manual']
    c = curl.Curl(url, silent = False, large = True)
    f = c.result
    result = []

    for l in f:

        l = l.strip().split('\t')

        if len(l) == 4:

            l[3] = ';'.join(
                map(
                    lambda x: x.split('-')[1].strip(),
                    filter(
                        lambda x: '-' in x,
                        l[3].replace('–', '-').split(',')
                    )
                )
            )

        l[0] = l[0].split('-')[0]
        l[1] = l[1].split('-')[0]
        result.append(
            NegatomeInteraction(*l)
        )

    return result
