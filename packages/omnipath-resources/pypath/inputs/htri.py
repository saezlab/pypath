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


def htri_interactions():

    HTRIInteraction = collections.namedtuple(
        'HTRIInteraction',
        [
            'entrez_tf',
            'genesymbol_tf',
            'entrez_target',
            'genesymbol_target',
            'pubmed',
        ]
    )

    c = curl.Curl(
        urls.urls['htri']['url'],
        init_url = urls.urls['htri']['init_url'],
        silent = False,
        follow = False,
        large = True,
    )

    data = c.result
    _ = next(c.result)

    return [
        HTRIInteraction(
            entrez_tf = fields[1],
            genesymbol_tf = fields[2],
            entrez_target = fields[3],
            genesymbol_target = fields[4],
            pubmed = fields[6],
        )
        for fields in
        (
            row.split(';') for row in data if row.strip()
        )
    ]
