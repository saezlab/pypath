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

import re
import collections

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.inputs.common as inputs_common


def macrophage_interactions():

    MacrophageInteraction = collections.namedtuple(
        'MacrophageInteraction',
        (
            'source_genesymbol',
            'target_genesymbol',
            'mechanism',
            'directed',
            'location',
            'pmid',
        )
    )


    url = urls.urls['macrophage']['url']
    c = curl.Curl(url, silent = False, large = True)
    fname = c.fileobj.name
    del c
    tbl = inputs_common.read_xls(fname)[5:]
    types = ['Protein', 'Complex']
    result = []

    for l in tbl:

        empty = {'', '-'}

        if len(l) > 11:

            if l[3].strip() in types and l[7].strip() in types:

                alist = _trim_gname(l[1])
                blist = _trim_gname(l[5])

                if len(alist) > 0 and len(blist) > 0:

                    for i in alist:

                        for j in blist:

                            if i != j and j not in empty and i not in empty:

                                pm = l[11].replace(',', '').strip().split('.')
                                pm = pm[0]

                                if not pm.startswith('INF'):

                                    directed = (
                                        '0'
                                            if l[9].strip() == 'Binding' else
                                        '1'
                                    )

                                    result.append(
                                        MacrophageInteraction(
                                            i,
                                            j,
                                            l[9].strip(),
                                            directed,
                                            l[10].strip(),
                                            pm,
                                        )
                                    )

    return result


def _trim_gname(gname):

    gname = re.sub(r'\[.*\]', '', re.sub(r'\(.*\)', '', gname))
    gname = re.sub(r'[A-Z]{0,1}[a-z]{1,}', '', gname)
    gname = gname.split(':')

    for i, g in enumerate(gname):

        gname[i] = gname[i].strip()

    return gname