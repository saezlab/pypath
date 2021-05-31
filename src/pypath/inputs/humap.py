#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2021
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

import itertools

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.utils.mapping as mapping
import pypath.internals.intera as intera


def humap_complexes():

    url = urls.urls['proteincomplexes']['url']
    c = curl.Curl(url, large = True)

    complexes = {}

    for l in c.result:

        l = l.strip().split()

        for uniprots in itertools.product(*(
            mapping.map_name(entrez, 'entrez', 'uniprot') for entrez in l
        )):

            cplex = intera.Complex(
                components = uniprots,
                sources = 'hu.MAP',
            )

            complexes[cplex.__str__()] = cplex

    return complexes
