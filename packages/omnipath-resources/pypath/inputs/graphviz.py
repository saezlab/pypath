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

import pypath.resources.urls as urls
import pypath.share.curl as curl


def graphviz_attrs():
    """
    Downloads the list of graphviz attributes from graphviz.org.
    Returns 3 dicts of dicts: graph_attrs, vertex_attrs and edge_attrs.
    """

    url = urls.urls['graphviz']['url']
    c = curl.Curl(url)
    html = c.result
    soup = bs4.BeautifulSoup(html, 'lxml')
    vertex_attrs = {}
    edge_attrs = {}
    graph_attrs = {}

    for tbl in soup.find_all('table'):

        if tbl.find('tr').text.strip().startswith('Name'):

            for r in tbl.find_all('tr'):

                r = r.find_all('td')

                if r:

                    usedby = r[1].text
                    this_attr = {
                        'type': r[2].text.strip(),
                        'default': r[3].text.strip(),
                        'min': r[4].text.strip(),
                        'notes': r[5].text.strip()
                    }
                    attr_name = r[0].text.strip()

                    if 'N' in usedby:

                        vertex_attrs[attr_name] = this_attr

                    if 'E' in usedby:

                        edge_attrs[attr_name] = this_attr

                    if 'G' in usedby:

                        graph_attrs[attr_name] = this_attr
            break

    return graph_attrs, vertex_attrs, edge_attrs
