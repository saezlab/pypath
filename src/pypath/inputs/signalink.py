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


def signalink_interactions(organism = 9606, exclude_secondary = True):
    """
    Reads and processes SignaLink3 interactions from local file.
    Returns list of interactions.
    """

    SignalinkInteraction = collections.namedtuple(
        'SignalinkInteraction',
        [
            'id_a',
            'id_b',
            'is_direct',
            'is_directed',
            'is_stimulation',
            'is_inhibition',
            'pathways_a',
            'pathways_b',
            'functions_a',
            'functions_b',
            'references',
            'resources',
        ]
    )

    signalink_sources = {
        'SLKv2.0',
        'SLKv2.1',
        'SLKv3.0',
        'SignaFish',
        'TCRcuration',
    }
    repref = re.compile(r'(?:.*:)?((?:[\w]+[^\s])?)\s?')
    nodes_pathways = {}
    nodes_organism = {}
    nodes_functions = {}
    interactions = []


    def get_value(field):

        return repref.sub('\\1', field)


    def get_values(field):

        return [get_value(f) for f in field.split('|')]


    url_nodes = urls.urls['signalink']['nodes']
    c_nodes = curl.Curl(url_nodes, silent = False, large = True)
    url_edges = urls.urls['signalink']['edges']
    c_edges = curl.Curl(url_edges, silent = False, large = True)

    _ = next(c_nodes.result)

    for l in c_nodes.result:

        l = l.strip('\n\r')

        if l:

            l = l.split('\t')
            _id = int(l[0])
            uniprot = get_value(l[1])
            pathways = l[4].split('|') if l[4] else []
            _organism = int(get_value(l[3]))
            nodes_pathways[uniprot] = pathways
            nodes_organism[uniprot] = _organism
            nodes_functions[uniprot] = l[6].split('|') if l[6] else []

    _ = next(c_edges.result)

    for l in c_edges.result:

        l = l.strip().split('\t')

        if exclude_secondary and not set(l[6].split('|')) & signalink_sources:

            continue

        id_a = get_value(l[0])
        id_b = get_value(l[1])

        if (
            nodes_organism[id_a] != organism or
            nodes_organism[id_b] != organism
        ):
            continue

        resources = [
            res for res in get_values(l[6])
            if not res.startswith('SLK')
        ]

        interaction_attrs = {
            tuple(iattr.split(':'))
            for iattr in l[5].split('|')
        }

        interactions.append(
            SignalinkInteraction(
                id_a = id_a,
                id_b = id_a,
                is_direct = ('is_direct', 'true') in interaction_attrs,
                is_directed = ('is_directed', 'true') in interaction_attrs,
                is_stimulation = (
                    ('MI', '0624(stimulation)') in interaction_attrs
                ),
                is_inhibition = (
                    ('MI', '0623(inhibition)') in interaction_attrs
                ),
                pathways_a = nodes_pathways[id_a],
                pathways_b = nodes_pathways[id_b],
                functions_a = nodes_functions[id_a],
                functions_b = nodes_functions[id_b],
                references = get_values(l[4]),
                resources = resources,
            )
        )

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
