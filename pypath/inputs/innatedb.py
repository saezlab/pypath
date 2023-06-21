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
import pypath.share.settings as settings


def innatedb_interactions(organism = 9606):

    InnatedbInteraction = collections.namedtuple(
        'InnatedbInteraction',
        (
            'source_uniprot',
            'source_genesymbol',
            'target_uniprot',
            'target_genesymbol',
            'pmid',
        )
    )

    url = urls.urls['innatedb']['url']
    headers = [settings.get('user_agent')]
    c = curl.Curl(url, silent = False, large = True, req_headers = headers)
    f = c.result
    result = []
    lnum = 0
    _ = next(c.result)

    for l in f:

        l = l.replace('\n', '').replace('\r', '')
        l = l.split('\t')
        specA = 0 if l[9] == '-' else int(l[9].split(':')[1].split('(')[0])
        specB = 0 if l[10] == '-' else int(l[10].split(':')[1].split('(')[0])

        if organism is None or (specA == organism and specB == organism):

            pm = l[8].replace('pubmed:', '')
            l = [l[4], l[5]]
            interaction = ()

            for ll in l:

                ll = ll.split('|')
                hgnc = ''
                uniprot = ''

                for lll in ll:

                    nm = lll.split(':')

                    if nm[0] == 'hgnc':

                        hgnc = nm[1].split('(')[0]

                    if nm[0] == 'uniprotkb' and len(nm[1]) == 6:

                        uniprot = nm[1]

                interaction += (uniprot, hgnc)

            interaction += (pm, )
            result.append(InnatedbInteraction(*interaction))

        lnum += 1

    f.close()

    return result
