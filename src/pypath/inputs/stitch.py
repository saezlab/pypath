#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2021
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: Dénes Türei (turei.denes@gmail.com)
#           Nicolàs Palacio
#           Olga Ivanova
#           Sebastian Lobentanzer
#           Ahmet Rifaioglu
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

import re
import collections

import pypath.resources.urls as urls
import pypath.share.curl as curl


def stitch_interactions(threshold = None):

    StitchInteraction = collections.namedtuple(
        'StitchInteraction',
        (
            'partner_a',
            'partner_b',
            'mechanism',
            'action',
            'score',
        ),
    )

    url = urls.urls['stitch']['actions']

    c = curl.Curl(url, silent = False, large = True)

    _ = next(c.result)

    sep = re.compile(r'[sm\.]')

    for l in c.result:

        if hasattr(l, 'decode'):

            l = l.decode('utf-8')

        l = l.strip().split('\t')

        score = int(l[5])

        if threshold is not None and score < threshold:
            continue

        try:
            a = sep.split(l[0])[1]
            b = sep.split(l[1])[1]

        except IndexError:
            print(l[1])

        if l[4] == 'f':
            a, b = b, a

        yield StitchInteraction(
            partner_a = a,
            partner_b = b,
            mechanism = l[2],
            action = l[3],
            score = int(l[5]),
        )
