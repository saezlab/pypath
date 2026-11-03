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


def _netbiol_interactions(database):

    record = collections.namedtuple(
        '%sInteraction' % database.capitalize(),
        (
            'source_uniprot',
            'target_uniprot',
            'is_direct',
            'is_directed',
            'effect',
        ) +
        (
            (
                'source_autophagy',
                'target_autophagy',
            )
            if database == 'arn' else
            ()
        ) +
        (
            'references',
        )
    )

    url = urls.urls[database]['url']
    c = curl.Curl(url, silent = True, large = False)

    return [
        record(*row.split(',')[:-1])
        for row in c.result.split('\n')
        if row
    ]


def arn_interactions():
    """
    Literature curated interactions between autophagy proteins from
    http://autophagyregulation.org/.
    """

    return _netbiol_interactions(database = 'arn')


def nrf2ome_interactions():
    """
    Literature curated interactions of the NRF2 pathway from
    http://nrf2ome.org/.
    """

    return _netbiol_interactions(database = 'nrf2ome')
