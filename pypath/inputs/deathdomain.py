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

import bs4

import pypath.resources.urls as urls
import pypath.share.curl as curl


def deathdomain_interactions():
    """
    Downloads HTML tables from the DeathDomain webpage and extracts
    the interactions.
    """

    result = []
    families = ['CARD', 'DD', 'DED', 'PYD']

    for fam in families:

        url = urls.urls['death']['url_dead'] % fam
        c = curl.Curl(url, silent = False)
        html = c.result

        soup = bs4.BeautifulSoup(html, 'lxml')

        d = {}
        for tab in soup.find_all('table', {'class': 'tab'}):
            for r in tab.find_all('tr'):
                cs = r.find_all('td')
                if len(cs) > 0:
                    i = {
                        'family': cs[0].find('a').text,
                        'A': cs[1].find('a').text,
                        'B': cs[3].find('a').text,
                        'met': cs[4].text if cs[4].text is not None else '',
                        'ref': cs[-1].find('a').text
                    }
                    if i['A'] not in d:
                        d[i['A']] = {}
                    if i['B'] not in d[i['A']]:
                        d[i['A']][i['B']] = {}
                    d[i['A']][i['B']]['family'] = i['family']
                    if 'met' not in d[i['A']][i['B']]:
                        d[i['A']][i['B']]['met'] = []
                    d[i['A']][i['B']]['met'].append(i['met'])
                    if 'ref' not in d[i['A']][i['B']]:
                        d[i['A']][i['B']]['ref'] = []
                    d[i['A']][i['B']]['ref'].append(i['ref'])

        for p1, v1 in iteritems(d):
            for p2, v2 in iteritems(v1):
                if p1 != p2:
                    result.append([
                        p1, p2, ';'.join(common.unique_list(v2['met'])),
                        ';'.join(common.unique_list(v2['ref']))
                    ])

    return result


def deathdomain_interactions_rescued():
    """
    Loads the DeathDomain interactions from rescued data.
    """

    url = urls.urls['death']['url_alive']

    c = curl.Curl(url, silent = False, large = True)

    _ = next(c.result)

    return [
        [i.strip() for i in line.split('\t')]
        for line in c.result
    ]
