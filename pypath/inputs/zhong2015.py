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

import collections

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.utils.mapping as mapping


def zhong2015_annotations():
    """
    From 10.1111/nyas.12776 (PMID 25988664).
    """

    types = {
        'i': 'iCAM',
        'm': 'matrix adhesion',
        'ag': 'axonal guidance',
        'aj': 'adherens junction',
        'c': 'cell-cell adhesion',
        'fa': 'focal adhesion',
        'tj': 'tight junction',
        'my': 'myelin interactions',
    }

    Zhong2015Annotation = collections.namedtuple(
        'Zhong2015Annotation',
        ['type'],
    )
    result = collections.defaultdict(set)

    url = urls.urls['zhong2015']['url']
    c = curl.Curl(url, silent = False, large = True)

    _ = next(c.result)

    for rec in c.result:
        rec = rec.strip().split('\t')

        uniprot = mapping.map_name0(rec[0], 'genesymbol', 'uniprot')

        if uniprot:
            result[uniprot].add(
                Zhong2015Annotation(type = types[rec[2]])
            )

    return dict(result)
