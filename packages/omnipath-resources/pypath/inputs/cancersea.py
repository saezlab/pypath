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

import os
import collections

import bs4

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.utils.mapping as mapping


def cancersea_annotations():
    """
    Retrieves genes annotated with cancer funcitonal states from the
    CancerSEA database.
    """

    CancerseaAnnotation = collections.namedtuple(
        'CancerseaAnnotation',
        [
            'state',
        ],
    )

    annotations = collections.defaultdict(set)

    url = urls.urls['cancersea']['rescued']
    c = curl.Curl(url, silent = False, large = False)
    soup = bs4.BeautifulSoup(c.result, 'html.parser')

    for row in soup.find_all('tbody')[1].find_all('tr'):

        state = row.find_all('td')[0].text
        url_end = row.find_all('td')[-1].find('a').attrs['href']
        url_end = url_end.rsplit('/', maxsplit = 1)[-1]
        data_url = urls.urls['cancersea']['rescued_data'] % url_end
        c = curl.Curl(data_url, silent = False, large = True)

        _ = next(c.result)

        for line in c.result:

            line = line.strip().split('\t')

            uniprots = mapping.map_name(line[1], 'genesymbol', 'uniprot')

            for uniprot in uniprots:

                annotations[uniprot].add(
                    CancerseaAnnotation(state = state)
                )

    return dict(annotations)
