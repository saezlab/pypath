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

import itertools

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.utils.mapping as mapping
import pypath.internals.intera as intera


def humap_complexes():

    url = urls.urls['humap']['humap_url']
    c = curl.Curl(url, large = True)

    complexes = {}

    for l in c.result:

        l = l.strip().split()

        for uniprots in itertools.product(*(
            mapping.map_name(entrez, 'entrez', 'uniprot')
            for entrez in l
        )):

            cplex = intera.Complex(
                components = uniprots,
                sources = 'hu.MAP',
            )

            complexes[cplex.__str__()] = cplex

    return complexes


def humap2_complexes(min_confidence = 0):

    url = urls.urls['humap']['humap2_url']
    c = curl.Curl(url, large = True)

    complexes = {}

    _ = next(c.result)

    for l in c.result:

        l = l.strip().split(',')

        confidence = int(l[1])

        if confidence < min_confidence:

            continue

        for uniprots in itertools.product(*(
            mapping.map_name(uniprot, 'uniprot', 'uniprot')
            for uniprot in l[2].split()
        )):

            if not uniprots:

                continue

            cplex = intera.Complex(
                components = uniprots,
                sources = 'hu.MAP2',
                attrs = {'humap2_confidence': confidence},
            )

            complexes[cplex.__str__()] = cplex

    return complexes