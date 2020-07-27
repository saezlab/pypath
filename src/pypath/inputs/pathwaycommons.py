#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2020
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#                  Nicolàs Palacio
#                  Olga Ivanova
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

# mapping input methods

from future.utils import iteritems

import collections

import pypath.share.curl as curl
import pypath.share.progress as progress
import pypath.resources.urls as urls
import pypath.share.common as common


def pathwaycommons_interactions(
        sources = None,
        types = None,
        sources_separated = True,
    ):

    result = {}
    interactions = []

    types = common.to_set(types)

    source_names = {
        'wp': 'WikiPathways',
        'kegg': 'KEGG',
        'bind': 'BIND',
        'intact': 'IntAct',
        'intact_complex': 'IntAct',
        'panther': 'PANTHER',
        'pid': 'NCI-PID',
        'reactome': 'Reactome',
        'dip': 'DIP',
        'hprd': 'HPRD',
        'inoh': 'INOH',
        'netpath': 'NetPath',
        'biogrid': 'BioGRID',
        'corum': 'CORUM',
        'psp': 'PhosphoSite',
    }

    directed = {
        'state-change',
        'controls-state-change-of',
        'controls-transport-of',
        'controls-phosphorylation-of',
    }

    sources = list(source_names.keys()) if sources is None else sources

    prg = progress.Progress(
        len(sources),
        'Processing PathwayCommons',
        1,
        percent = False,
    )

    url = urls.urls['pwcommons']['url']

    for s in sources:

        prg.step()
        surl = url % s
        c = curl.Curl(surl, silent = False, large = True)

        for l in c.result:

            if hasattr(l, 'decode'):

                l = l.decode('ascii')

            l = l.strip().split('\t')

            if types is None or l[1] in types:

                if sources_separated:

                    l.append(source_names[s])
                    interactions.append(l)

                else:

                    pair = (l[0], l[2])

                    if pair not in result:

                        result[pair] = [set(), set(), 0]

                    result[pair][0].add(source_names[s])
                    result[pair][1].add(l[1])

                    if l[1] in directed:

                        result[pair][2] = 1

    if not sources_separated:
        for pair, details in iteritems(result):
            interactions.append([
                pair[0], pair[1], ';'.join(details[0]), ';'.join(details[1]),
                str(details[2])
            ])

    return interactions