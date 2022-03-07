#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2022
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: Dénes Türei (turei.denes@gmail.com)
#           Nicolàs Palacio
#           Olga Ivanova
#           Sebastian Lobentanzer
#           Ahmet Rifaioglu
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

import collections

import pypath.resources.urls as urls
import pypath.share.curl as curl


def biogrid_interactions(organism = 9606, htp_limit = 1, ltp = True):
    """
    Downloads and processes Physical multi-validated BioGRID interactions.
    Keeps only the "low throughput" interactions.
    Returns list of interactions.

    @organism : int
        NCBI Taxonomy ID of organism.
    @htp_limit : int
        Exclude interactions only from references
        cited at more than this number of interactions.
    """

    BiogridInteraction = collections.namedtuple(
        'BiogridInteraction',
        (
            'partner_a',
            'partner_b',
            'pmid',
        ),
    )

    organism = str(organism)
    interactions = []
    refc = []
    url = urls.urls['biogrid']['mv']
    c = curl.Curl(url, silent = False, large = True)
    f = next(iter(c.result.values()))
    nul = f.readline()

    for l in f:

        l = l.split('\t')

        if len(l) > 17:

            if (
                (
                    l[17].startswith('Low') or
                    not ltp
                ) and (
                    l[15] == organism and
                    l[16] == organism
                )
            ):

                interactions.append(
                    BiogridInteraction(
                        partner_a = l[7],
                        partner_b = l[8],
                        pmid = l[14]
                    )
                )
                refc.append(l[14])

    refc = collections.Counter(refc)

    if htp_limit is not None:

        interactions = [i for i in interactions if refc[i[2]] <= htp_limit]

    return interactions


def biogrid_all_interactions(organism = 9606, htp_limit = 1, ltp = True):

    """
    Downloads and processes all BioGRID interactions.
    Keeps only the "low throughput" interactions.
    Returns list of interactions.

    @organism : int
        NCBI Taxonomy ID of organism.
    @htp_limit : int
        Exclude interactions only from references
        cited at more than this number of interactions.
    """


    BiogridInteraction = collections.namedtuple(
        'BiogridInteraction',
        (
            'partner_a',
            'partner_b',
            'pmid',
            'experimental_system',
            'experimental_system_type',
            'throughput',
            'htp_score',
            'multi_validated',
        ),
    )

    organism = str(organism)
    interactions = []
    refc = []
    mv_dict=collections.defaultdict(list)
    for i in biogrid_interactions(organism, htp_limit, ltp):
        mv_dict[i.partner_a].append(i.partner_b)
    url = urls.urls['biogrid']['all']
    c = curl.Curl(url, silent = False, large = True)
    f = next(iter(c.result.values()))
    nul = f.readline()

    for l in f:

        l = l.split('\t')

        if len(l) > 17:

            if (
                (
                    l[17].startswith('Low') or
                    not ltp
                ) and (
                    l[15] == organism and
                    l[16] == organism
                )
            ):
                if l[8] in mv_dict[l[7]] or l[7] in mv_dict[l[8]]:
                    mv="+"
                else:
                    mv=""
                interactions.append(
                    BiogridInteraction(
                        partner_a = l[7],
                        partner_b = l[8],
                        pmid = l[14],
                        experimental_system= l[11],
                        experimental_system_type= l[12],
                        throughput= l[17],
                        htp_score= l[18],
                        multi_validated= mv
                    )
                )
                refc.append(l[14])

    refc = collections.Counter(refc)

    if htp_limit is not None:

        interactions = [i for i in interactions if refc[i[2]] <= htp_limit]

    return interactions
