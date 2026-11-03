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


def encode_tf_mirna_interactions():

    EncodeInteraction = collections.namedtuple(
        'EncodeInteraction',
        (
            'tf_genesymbol',
            'mirna',
        ),
    )

    url = urls.urls['encode']['tf-mirna']
    c = curl.Curl(
        url,
        silent = False,
        large = True,
        encoding = 'ascii',
    )

    result = []

    for l in c.result:

        l = l.strip().split()

        if l[1] == '(TF-miRNA)':

            result.append(
                EncodeInteraction(l[0], l[2])
            )

    return result
