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
import collections

import pypath.resources.urls as urls
import pypath.share.curl as curl


def ccmap_interactions(organism = 9606):
    """
    Downloads and processes CancerCellMap.
    Returns list of interactions.

    @organism : int
        NCBI Taxonomy ID to match column #7 in nodes file.
    """

    CancercellmapInteraction = collections.namedtuple(
        'CancercellmapInteraction',
        (
            'source_uniprot',
            'target_uniprot',
            'directed',
            'references',
        ),
    )

    organism = '%u' % organism
    interactions = []
    nodes_url = urls.urls['ccmap']['nodes']
    edges_url = urls.urls['ccmap']['edges']

    c = curl.Curl(
        nodes_url,
        silent = False,
        files_needed = ['cell-map-node-attributes.txt'],
    )
    nodes = c.result

    c = curl.Curl(
        edges_url,
        silent = False,
        files_needed = ['cell-map-edge-attributes.txt'],
    )
    edges = c.result

    nodes = dict(
        map(
            lambda l: (l[1], l[2].split(':')),
            filter(
                lambda l: l[5] == 'protein' and l[6] == organism,
                filter(
                    lambda l: len(l) == 7,
                    map(
                        lambda l: l.strip().split('\t'),
                        nodes['cell-map-node-attributes.txt'].split('\n')[1:]
                    )
                )
            )
        )
    )

    edges = filter(
        lambda l: len(l) == 7,
        map(
            lambda l: l.strip().split('\t'),
            edges['cell-map-edge-attributes.txt'].split('\n')[1:]
        )
    )

    for e in edges:

        if e[1] != 'IN_SAME_COMPONENT' and e[3] in nodes and e[4] in nodes:

            for src, tgt in itertools.product(nodes[e[3]], nodes[e[4]]):

                interactions.append(
                    CancercellmapInteraction(
                        source_uniprot = src,
                        target_uniprot = tgt,
                        directed = e[1] == 'STATE_CHANGE',
                        references = e[6].strip(';').replace('PUBMED:', ''),
                    )
                )

    return interactions
