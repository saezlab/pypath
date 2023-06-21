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

import bs4

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.utils.mapping as mapping


def get_integrins():
    """
    Returns a set of the UniProt IDs of the human integrins from
    Table 1 of Takada et al 2007 (10.1186/gb-2007-8-5-215).
    """

    url = urls.urls['integrins']['url']

    req_headers = [
        'Host: www.ncbi.nlm.nih.gov',
        'User-Agent: Mozilla/5.0 (X11; Linux x86_64; rv:67.0) '\
            'Gecko/20100101 Firefox/67.0',
        'Accept: text/html,application/xhtml+xml,'
            'application/xml;q=0.9,*/*;q=0.8',
        'Accept-Language: en-US,en;q=0.5',
        'Connection: keep-alive',
        'Upgrade-Insecure-Requests: 1',
        'Pragma: no-cache',
        'Cache-Control: no-cache',
    ]

    c = curl.Curl(
        url, silent = False, req_headers = req_headers, large = True,
    )
    soup = bs4.BeautifulSoup(c.fileobj.read(), 'lxml')

    integrins = []

    rows = soup.find_all('tr')

    for tr in rows[1:]:
        cells = [td for td in tr.find_all('td')]
        integrins.append(cells[-1].text.split('}')[-1])

    return mapping.map_names(integrins, 'uniprot', 'uniprot')