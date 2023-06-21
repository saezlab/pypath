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

from future.utils import iteritems

import re
import bs4

import pypath.utils.taxonomy as taxonomy
import pypath.share.curl as curl
import pypath.resources.urls as urls


def phosphoelm_enzyme_substrate(organism = 9606, ltp_only = True):
    """
    Downloads kinase-substrate interactions from phosphoELM.
    Returns list of dicts.

    :param int organism: NCBI Taxonomy ID.
    :param bool ltp_only: Include only low-throughput interactions.
    """

    result = []
    non_digit = re.compile(r'[^\d.-]+')

    if organism is None:
        _organism = None
    elif organism in taxonomy.phosphoelm_taxids:
        _organism = taxonomy.phosphoelm_taxids[organism]
    else:
        sys.stdout.write('\t:: Unknown organism: `%u`.\n' % organism)
        return []

    url = urls.urls['p_elm']['url']
    c = curl.Curl(url, silent = False)
    data = c.result
    data = [
        n for d, n in iteritems(data)
        if d.startswith(urls.urls['p_elm']['psites'])
    ]
    data = data[0] if len(data) > 0 else ''
    data = [l.split('\t') for l in data.split('\n')]
    kinases = phosphoelm_kinases()
    del data[0]

    for l in data:

        if (
            len(l) == 9 and (
                l[7] == _organism or
                _organism is None
            ) and (
                not ltp_only or
                l[6] == 'LTP'
            )
        ):

            l[1] = 1 if '-' not in l[0] else int(l[0].split('-')[1])
            l[0] = l[0].split('-')[0]
            del l[-1]

            if len(l[5]) > 0 and l[5] in kinases:
                kinase = kinases[l[5]]

                result.append({
                    'instance': None,
                    'isoform': l[1],
                    'resaa': l[3],
                    'resnum': int(non_digit.sub('', l[2])),
                    'start': None,
                    'end': None,
                    'substrate': l[0],
                    'kinase': kinase,
                    'references': l[4].split(';'),
                    'experiment': l[6],
                    'organism': l[7]
                })

    return result


def phosphoelm_interactions(organism = 'Homo sapiens'):

    result = []
    data = phosphoelm_enzyme_substrate(ltp_only = True)

    for l in data:

        result.append([
            l['kinase'],
            l['substrate'],
            ';'.join(l['references']),
            l['organism'],
        ])

    return result


def phosphoelm_kinases():

    result = {}
    url = urls.urls['p_elm_kin']['url']
    c = curl.Curl(url, silent = False)
    data = c.result
    soup = bs4.BeautifulSoup(data, 'html.parser')

    for row in soup.find_all('table')[1].find_all('tr'):

        thisRow = [x.text for x in row.find_all('td')]

        if len(thisRow) > 2 and len(thisRow[2].strip()) > 0:

            result[thisRow[0]] = thisRow[2].strip()

    return result
