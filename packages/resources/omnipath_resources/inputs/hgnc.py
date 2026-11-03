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

def hgnc_genegroups():

    HGNCGeneGroupAnnotation = collections.namedtuple(
        'HGNCGeneGroupAnnotation',
        ['mainclass'],
    )


    result = collections.defaultdict(set)

    url = urls.urls['hgnc']['groups']
    c = curl.Curl(url, large = True, silent = False)

    _ = next(c.result)

    for rec in c.result:

        rec = rec.split('\t')
        uniprots = {u.strip() for u in rec[2].split(',')}
        uniprots.discard('')

        if not uniprots:
            continue

        uniprots = mapping.map_names(uniprots, 'uniprot', 'uniprot')

        if not uniprots:
            continue

        groups = rec[3].split('|')

        for group in groups:
            group = group.strip()

            if group:
                for uniprot in uniprots:
                    result[uniprot].add(
                        HGNCGeneGroupAnnotation(mainclass = group)
                    )

    return dict(result)
