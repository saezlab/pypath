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

from typing import List, Optional
from numbers import Number

import collections

import pypath.resources.urls as urls
import pypath.share.curl as curl


def biogrid_interactions(
        organism: int = 9606,
        htp_limit: Optional[Number] = 1,
        ltp: bool = True,
    ) -> List[tuple]:
    """
    Downloads and processes Physical multi-validated BioGRID interactions.
    Keeps only the "low throughput" interactions.
    Returns list of interactions.

    Args
        organism: NCBI Taxonomy ID of organism.
        htp_limit: Exclude interactions only from references cited at more
        than this number of interactions.
    """

    BiogridPhysicalInteraction = collections.namedtuple(
        'BiogridPhysicalInteraction',
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
    c = curl.Curl(url, silent = False, large = True, slow = True)
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
                    organism == 'None' or (
                        l[15] == organism and
                        l[16] == organism
                    )
                )
            ):

                interactions.append(
                    BiogridPhysicalInteraction(
                        partner_a = l[7],
                        partner_b = l[8],
                        pmid = l[14],
                    )
                )
                refc.append(l[14])

    refc = collections.Counter(refc)

    if htp_limit is not None:

        interactions = [i for i in interactions if refc[i[2]] <= htp_limit]

    return interactions


def biogrid_all_interactions(
        organism: int = 9606,
        htp_limit: Optional[Number] = 1,
        ltp: bool = True,
    ) -> List[tuple]:
    """
    Downloads and processes all BioGRID interactions.
    Keeps only the "low throughput" interactions.
    Returns list of interactions.

    Args
        organism: NCBI Taxonomy ID of organism.
        htp_limit: Exclude interactions only from references cited at
            more than this number of interactions.
    """


    BiogridInteraction = collections.namedtuple(
        'BiogridInteraction',
        (
            'partner_a',
            'partner_b',
            'pmid',
            'experimental_system',
            'experimental_system_type',
            'ltp',
            'htp_score',
            'multi_validated',
            'tax_a',
            'tax_b',
        ),
    )

    organism = str(organism)
    interactions = []
    refc = []
    mv_dict=collections.defaultdict(list)

    for i in biogrid_interactions(organism, htp_limit, ltp):

        mv_dict[i.partner_a].append(i.partner_b)

    url = urls.urls['biogrid']['all']
    c = curl.Curl(url, silent = False, large = True, slow = True)
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
                    organism == 'None' or (
                        l[15] == organism and
                        l[16] == organism
                    )
                )
            ):

                mv = l[8] in mv_dict[l[7]] or l[7] in mv_dict[l[8]]

                interactions.append(
                    BiogridInteraction(
                        partner_a = l[7],
                        partner_b = l[8],
                        pmid = l[14],
                        experimental_system = l[11],
                        experimental_system_type = l[12],
                        ltp = 'Low T' in l[17],
                        htp_score = None if l[18] == '-' else float(l[18]),
                        multi_validated = mv,
                        tax_a = l[15],
                        tax_b = l[16],
                    )
                )

                refc.append(l[14])

    refc = collections.Counter(refc)

    if htp_limit is not None:

        interactions = [i for i in interactions if refc[i[2]] <= htp_limit]

    return interactions
