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


def humancellmap_annotations():

    HumancellmapAnnotation = collections.namedtuple(
        'HumancellmapAnnotation',
        [
            'localization',
            'method',
        ],
    )

    url = urls.urls['humancellmap']['preys_url']

    c = curl.Curl(url, large = True, silent = False)

    result = collections.defaultdict(set)

    _ = next(c.result)

    for l in c.result:

        l = l.strip().split('\t')

        if len(l) < 10:

            continue

        for method, i_loc in (('NMF', 2), ('SAFE', 4)):

            for localization in l[i_loc].split(','):

                localization = localization.strip()

                if localization in {'-', 'no prediction'}:

                    continue

                for uniprot in mapping.map_name(l[9], 'uniprot', 'uniprot'):

                    result[uniprot].add(
                        HumancellmapAnnotation(
                            localization = localization,
                            method = method,
                        )
                    )

    return dict(result)

