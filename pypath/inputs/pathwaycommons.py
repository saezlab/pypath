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

import collections

import pypath.share.curl as curl
import pypath.share.progress as progress
import pypath.resources.urls as urls
import pypath.share.common as common


PathwayCommonsInteraction = collections.namedtuple(
    'PathwayCommonsInteraction',
    [
        'id_a',
        'interaction_type',
        'id_b',
        'resource',
    ],
)


PathwayCommonsResource = collections.namedtuple(
    'PathwayCommonsResource',
    [
        'name',
        'pc_label',
        'version',
    ]
)
PathwayCommonsResource.__new__.__defaults__ = (12,)


pathwaycommons_resources = {
    PathwayCommonsResource('WikiPathways', 'wp', 11),
    PathwayCommonsResource('KEGG', 'kegg'),
    PathwayCommonsResource('BIND', 'bind'),
    PathwayCommonsResource('IntAct', 'intact'),
    PathwayCommonsResource('IntAct', 'intact_complex'),
    PathwayCommonsResource('PANTHER', 'panther'),
    PathwayCommonsResource('NCI-PID', 'pid'),
    PathwayCommonsResource('Reactome', 'reactome'),
    PathwayCommonsResource('DIP', 'dip'),
    PathwayCommonsResource('HPRD', 'hprd'),
    PathwayCommonsResource('INOH', 'inoh'),
    PathwayCommonsResource('NetPath', 'netpath'),
    PathwayCommonsResource('BioGRID', 'biogrid'),
    PathwayCommonsResource('CORUM', 'corum'),
    PathwayCommonsResource('PhosphoSite', 'psp'),
}


pathwaycommons_directed_types = {
    'state-change',
    'controls-state-change-of',
    'controls-transport-of',
    'controls-phosphorylation-of',
}


def pathwaycommons_interactions(
        resources = None,
        types = None,
        by_interaction = False,
        version = 12,
    ):

    interactions = collections.defaultdict(set) if by_interaction else []

    types = common.to_set(types)

    resources = {
        res.lower()
        for res in (
            common.to_list(resources) or (
                pc_res.name
                for pc_res in pathwaycommons_resources
            )
        )
    }

    prg = progress.Progress(
        len(resources),
        'Processing PathwayCommons',
        1,
        percent = False,
    )

    url = urls.urls['pwcommons']['url']

    for resource in pathwaycommons_resources:

        if not resources & {resource.pc_label, resource.name.lower()}:

            continue

        prg.step()
        _version = min(resource.version, version)
        resource_url = url % (_version, _version, resource.pc_label)
        c = curl.Curl(resource_url, silent = False, large = True)

        for l in c.result:

            if hasattr(l, 'decode'):

                l = l.decode('ascii')

            l = l.strip('\n\r').split('\t')

            if not types or l[1] in types:

                if by_interaction:

                    a_b = (l[0], l[1], l[2])
                    b_a = (l[2], l[1], l[0])

                    directed = l[1] in pathwaycommons_directed_types

                    key = (
                        b_a
                            if (
                                a_b not in interactions and
                                not directed and
                                b_a in interactions
                            ) else
                        a_b
                    )

                    interactions[key].add(
                        PathwayCommonsInteraction(
                            *key,
                            resource = resource.name
                        )
                    )

                else:

                    l.append(resource.name)
                    interactions.append(PathwayCommonsInteraction(*l))

    return interactions


def _create_single_resource_method(resource):

    def _pc_single_resource(**kwargs):

        kwargs['resources'] = resource.name

        return pathwaycommons_interactions(**kwargs)


    return _pc_single_resource


for res in pathwaycommons_resources:

    method_name = 'pathwaycommons_%s_interactions' % res.name.lower()
    globals()[method_name] = _create_single_resource_method(res)