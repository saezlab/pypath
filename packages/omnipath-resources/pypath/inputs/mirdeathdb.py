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

import pypath.resources.urls as urls
import pypath.share.curl as curl


def mirdeathdb_interactions():
    """
    Retrieves literature curated miRNA-target gene interactions from
    miRDeathDB (https://www.nature.com/articles/cdd201287).
    """

    MirdeathdbInteraction = collections.namedtuple(
        'MirdeathdbInteraction',
        (
            'mirna',
            'target_entrez',
            'organism',
            'pubmed',
            'function',
        ),
    )

    url = urls.urls['mirdeathdb']['url_rescued']
    c = curl.Curl(url, silent = False, large = True)

    _ = next(c.result)

    for l in c.result:

        l = l.strip().split('\t')

        if len(l) < 11:

            continue

        mirnas = l[2].replace('"', '').split(',')
        organism = int(l[9])
        pubmed = l[8]
        geneid = l[10]
        function = '%s_%s' % (l[4], l[5])

        for mirna in mirnas:

            yield MirdeathdbInteraction(
                mirna = mirna.strip(),
                target_entrez = geneid,
                organism = organism,
                pubmed = pubmed,
                function = function,
            )
