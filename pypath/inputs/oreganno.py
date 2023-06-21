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

import re
import collections

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.utils.taxonomy as taxonomy


def oreganno_raw():
    """
    Downloads TF-target data from the ORegAnnO database.

    Yields:
        Tuples of raw records.
    """

    url = urls.urls['oreganno']['url']
    c = curl.Curl(url, silent = False, large = True, slow = True)
    data = c.result
    _ = next(data)

    for l in data:

        if not l:

            continue

        yield tuple(x.strip() for x in l.split('\t'))


def oreganno_interactions(organism = 9606):
    """
    Downloads TF-target interactions from the ORegAnnO database.

    Yields:
        Named tuples of TF, target and literature references.
    """

    OregannoInteraction = collections.namedtuple(
        'OregannoInteraction',
        ('tf', 'target', 'pmid'),
    )

    taxids = taxonomy.phosphoelm_taxids

    if organism in taxids:

        organism = taxids[organism]

    nsep = re.compile(r'([-A-Za-z0-9]{3,})[\s/\(]*.*')
    nrem = re.compile(r'[-/]')

    for l in oreganno_raw():

        if (l[1] == organism and
            l[3] == 'TRANSCRIPTION FACTOR BINDING SITE' and
            l[2] == 'POSITIVE OUTCOME' and
            l[4] != 'N/A' and
            l[7] != 'N/A'
        ):

            yield OregannoInteraction(
                tf = (
                    l[7]
                        if len(l[7]) < 3 else
                    nrem.sub('', nsep.findall(l[7])[0])
                ),
                target = (
                    l[4]
                        if len(l[4]) < 3 else
                    nrem.sub('', nsep.findall(l[4])[0])
                ),
                pmid = l[11] if l[11] != 'N/A' else '',
            )
