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

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.common as common


def cpdb_interactions(exclude = None):
    """
    Interactions from ConsensusPathDB.

    Args
        exclude (set): A set of resource names to exclude.
    """

    exclude = common.to_set(exclude)
    result = []
    url = urls.urls['cpdb']['url']
    c = curl.Curl(url, silent = False)
    data = c.result
    data = [
        x.split('\t')
        for x in data.split('\n')
        if not x.startswith('#') and len(x) > 0
    ]

    for l in data:

        participants = l[2].split(',')

        if len(participants) == 2:

            if (
                not exclude or
                set(l[0].split(',')) - exclude
            ):

                result.append([
                    participants[0],
                    participants[1],
                    l[0],
                    l[1],
                ])

    return result


def cpdb_interactions_ltp():
    """
    Low-throughput interactions from ConsensusPathDB. Calls
    ``cpdb_interactions`` with excluding HPRD, BioGRID, PhosphoPOINT,
    MINT, BIND and IntAct.

    Args
        exclude (set): A set of resource names to exclude.
    """

    return cpdb_interactions(
        exclude = {
            'HPRD',
            'BioGRID',
            'PhosphoPOINT',
            'MINT',
            'BIND',
            'IntAct',
        }
    )
