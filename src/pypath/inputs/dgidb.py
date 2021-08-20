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

import csv
import collections

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.utils.mapping as mapping


def dgidb_annotations():
    """
    Downloads druggable protein annotations from DGIdb.
    """

    DgidbAnnotation = collections.namedtuple(
        'DgidbAnnotation',
        ['category'],
    )


    url = urls.urls['dgidb']['categories']
    c = curl.Curl(url = url, silent = False, large = True)
    data = csv.DictReader(c.result, delimiter = '\t')

    result = collections.defaultdict(set)

    for rec in data:

        uniprots = mapping.map_name(
            rec['entrez_gene_symbol'],
            'genesymbol',
            'uniprot',
        )

        for uniprot in uniprots:
            result[uniprot].add(
                DgidbAnnotation(
                    category = rec['category']
                )
            )

    return result
