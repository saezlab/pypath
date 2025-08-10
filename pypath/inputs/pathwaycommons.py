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
        version = 14,  # Updated to use v14 by default
    ):

    interactions = collections.defaultdict(set) if by_interaction else []

    types = common.to_set(types)

    # For v14+, use unified file with source information; for older versions, use individual resource files
    if version >= 14:
        # Use unified txt file from v14 which has source information
        unified_url = (
            'https://download.baderlab.org/PathwayCommons/PC2/'
            f'v{version}/pc-hgnc.txt.gz'
        )
        c = curl.Curl(unified_url, silent = False, large = True, slow = True)

        if c.result is None:
            return interactions

        # Skip header line
        header_line = next(c.result, None)

        # Map resource names for filtering
        resource_map = {
            'wikipathways': 'WikiPathways',
            'kegg': 'KEGG',
            'biogrid': 'BioGRID',
            'reactome': 'Reactome',
            'intact': 'IntAct',
            'pid': 'NCI-PID',
            'hprd': 'HPRD',
            'dip': 'DIP',
            'bind': 'BIND',
            'corum': 'CORUM',
            'panther': 'PANTHER',
            'netpath': 'NetPath',
            'inoh': 'INOH',
        }

        # Convert requested resources to expected format
        if resources:
            filtered_resources = set()
            for res in resources:
                res_lower = res.lower()
                # Add both the mapped name and original name for matching
                if res_lower in resource_map:
                    filtered_resources.add(resource_map[res_lower])
                filtered_resources.add(res_lower)
                filtered_resources.add(res)
        else:
            filtered_resources = None

        for l in c.result:

            if hasattr(l, 'decode'):
                l = l.decode('ascii')

            l = l.strip('\n\r').split('\t')

            if len(l) >= 4:  # Now we have: [A, interaction_type, B, data_source, pubmed, pathways, ...]
                interaction_type = l[1]
                data_sources = l[3].split(';') if l[3] else []

                # Filter by data source if resources specified
                if filtered_resources:
                    source_match = any(
                        src.lower() in filtered_resources or src in filtered_resources
                        for src in data_sources
                    )
                    if not source_match:
                        continue

                # Filter by interaction type
                if types and interaction_type not in types:
                    continue

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
                            resource = ';'.join(data_sources)  # Use actual data sources
                        )
                    )

                else:

                    interactions.append(PathwayCommonsInteraction(
                        l[0], l[1], l[2], ';'.join(data_sources)
                    ))

    else:
        # Use old individual resource files for v12 and below
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

            if c.result is None:
                continue

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
