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

import collections

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.common as common


def string_effects(
    ncbi_tax_id = 9606,
    stimulation = ('activation',),
    inhibition = ('inhibition',),
    exclude = ('expression',),
    score_threshold = 0,
):

    StringInteraction = collections.namedtuple(
        'StringInteraction',
        (
            'source',
            'target',
            'effect',
        ),
    )

    effects = []

    stimulation = common.to_set(stimulation)
    inhibition = common.to_set(inhibition)
    exclude = common.to_set(exclude)

    url = urls.urls['string']['actions'] % ncbi_tax_id
    c = curl.Curl(url, silent = False, large = True)
    _ = next(c.result)

    for l in c.result:

        if hasattr(l, 'decode'):

            l = l.decode('ascii')

        l = l.strip().split('\t')

        if l and l[4] == 't' and int(l[6]) >= score_threshold:

            effect = (
                '+'
                    if l[2] in stimulation else
                '-'
                    if l[2] in inhibition else
                '*'
                    if l[2] not in exclude else
                None
            )

            source = l[0].split('.')[1] if l[5] == 't' else l[1].split('.')[1]
            target = l[1].split('.')[1] if l[5] == 't' else l[0].split('.')[1]

            if effect is not None:

                effects.append(
                    StringInteraction(
                        source = source,
                        target = target,
                        effect = effect,
                    )
                )

    return effects
