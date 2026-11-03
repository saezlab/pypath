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

from future.utils import iteritems

import sys
import re
import collections
import itertools

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.share.common as common
import pypath.utils.mapping as mapping


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
            'effect',
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
            pathways = (
                l[4].replace('/Wingless', '').split('|') if l[4] else []
            )
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

        effects = []

        if ('MI', '0624(stimulation)') in interaction_attrs:

            effects.append(1)

        if ('MI', '0623(inhibition)') in interaction_attrs:

            effects.append(-1)

        if not effects:

            effects.append(0)

        for effect in effects:

            interactions.append(
                SignalinkInteraction(
                    id_a = id_a,
                    id_b = id_b,
                    is_direct = ('is_direct', 'true') in interaction_attrs,
                    is_directed = ('is_directed', 'true') in interaction_attrs,
                    effect = effect,
                    pathways_a = nodes_pathways[id_a],
                    pathways_b = nodes_pathways[id_b],
                    functions_a = nodes_functions[id_a],
                    functions_b = nodes_functions[id_b],
                    references = get_values(l[4]),
                    resources = resources,
                )
            )

    return interactions


def signalink_annotations(organism = 9606):

    SignalinkPathway = collections.namedtuple(
        'SignalinkPathway',
        [
            'pathway',
        ]
    )

    SignalinkFunction = collections.namedtuple(
        'SignalinkFunction',
        [
            'function',
        ]
    )


    result = {
        'pathway': collections.defaultdict(set),
        'function': collections.defaultdict(set),
    }

    interactions = signalink_interactions(organism = organism)

    for i in interactions:

        for postfix in ('_a', '_b'):

            _id = getattr(i, 'id%s' % postfix)

            for uniprot in mapping.map_name(_id, 'uniprot', 'uniprot'):

                for attr, record in zip(
                    ('pathway', 'function'),
                    (SignalinkPathway, SignalinkFunction),
                ):

                    values = getattr(i, '%ss%s' % (attr, postfix))

                    for value in values:

                        result[attr][uniprot].add(
                            record(value)
                        )

    return dict((k, dict(v)) for k, v in iteritems(result))


def signalink_pathway_annotations(organism = 9606):

    return signalink_annotations(organism = organism)['pathway']


def signalink_function_annotations(organism = 9606):

    return signalink_annotations(organism = organism)['function']
