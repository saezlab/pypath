#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2020
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

from future.utils import iteritems

import sys
import re
import collections

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.share.common as common


def signalink_interactions():
    """
    Reads and processes SignaLink3 interactions from local file.
    Returns list of interactions.
    """

    repar = re.compile(r'.*\(([a-z\s]+)\)')
    repref = re.compile(r'(?:.*:)?((?:[\w]+[^\s])?)\s?')
    notNeeded = set(['acsn', 'reactome'])
    nodes = {}
    interactions = []

    def _get_attr(attrs, attrName):
        return _process_attr(attrs[attrName]) if attrName in attrs else ''

    def _process_attr(attr):
        m = repar.match(attr)

        if m is not None:
            return m.groups()[0]

        else:
            return attr

    url_nodes = urls.urls['signalink']['nodes']
    c_nodes = curl.Curl(url_nodes, silent = False, large = True)
    url_edges = urls.urls['signalink']['edges']
    c_edges = curl.Curl(url_edges, silent = False, large = True)

    for l in c_nodes.result:
        if len(l) > 0:
            l = l.split('\t')
            _id = int(l[0])
            uniprot = repref.sub('\\1', l[1])
            pathways = [
                pw.split(':')[-1].strip() for pw in l[4].split('|')
                if pw.split(':')[0] not in notNeeded
            ]
            nodes[_id] = [uniprot, pathways]

    lPrev = None

    for l in c_edges.result:
        l = l.strip().split('\t')

        if lPrev is not None:
            l = lPrev + l[1:]
            lPrev = None

        if len(l) == 13:
            if l[-1] == '0':

                dbs = [
                    _process_attr(db.split(':')[-1])
                    for db in l[9].replace('"', '').split('|')
                ]
                dbs = list(set(dbs) - notNeeded)

                if len(dbs) == 0:
                    continue

                idSrc = int(l[1])
                idTgt = int(l[2])

                uniprotSrc = repref.sub('\\1', l[3])
                uniprotTgt = repref.sub('\\1', l[4])

                if not uniprotSrc or not uniprotTgt:

                    continue

                refs = [ref.split(':')[-1] for ref in l[7].split('|')]
                attrs = dict(
                    tuple(attr.strip().split(':', 1))
                    for attr in l[8].replace('"', '').split('|'))
                interactions.append([
                    uniprotSrc, uniprotTgt, ';'.join(refs), ';'.join(dbs),
                    _get_attr(attrs, 'effect'),
                    _get_attr(attrs, 'is_direct'),
                    _get_attr(attrs, 'is_directed'),
                    _get_attr(attrs, 'molecular_background'),
                    ';'.join(nodes[idSrc][1]), ';'.join(nodes[idTgt][1])
                ])

        else:
            lPrev = l

    return interactions


def signalink_pathway_annotations():

    SignalinkPathway = collections.namedtuple(
        'SignalinkPathway',
        ['pathway', 'core'],
    )


    result = collections.defaultdict(set)

    interactions = signalink_interactions()

    for i in interactions:
        for idx in (0, 1):
            for pathway in i[idx + 8].split(';'):

                core = 'non-core' not in pathway
                pathway = (
                    pathway.split('(')[0].strip().replace('/Wingless', '')
                )

                for uniprot in mapping.map_name(i[idx], 'uniprot', 'uniprot'):

                    result[uniprot].add(
                        SignalinkPathway(pathway = pathway, core = core)
                    )

    return result
