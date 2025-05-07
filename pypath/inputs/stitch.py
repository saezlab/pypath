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

from typing import List, Literal, Union
from numbers import Number

import re
import collections

import pypath.resources.urls as urls
import pypath.share.curl as curl


def stitch_actions_interactions(threshold: Number = None) -> List[tuple]:

    StitchActionsInteraction = collections.namedtuple(
        'StitchActionsInteraction',
        (
            'partner_a',
            'partner_b',
            'mechanism',
            'action',
            'score',
        ),
    )

    url = urls.urls['stitch']['actions']

    c = curl.Curl(url, silent = False, large = True, slow = True)

    _ = next(c.result)  

    sep = re.compile(r'[sm\.]')

    for l in c.result:

        if hasattr(l, 'decode'):

            l = l.decode('utf-8')

        l = l.strip().split('\t')

        score = int(l[5])

        if threshold is not None and score < threshold:

            continue

        a = sep.split(l[0])[1]
        b = sep.split(l[1])[1]

        if l[4] == 'f':

            a, b = b, a

        yield StitchActionsInteraction(
            partner_a = a,
            partner_b = b,
            mechanism = l[2],
            action = l[3] or None,
            score = int(l[5]),
        )


def stitch_links_interactions(
        ncbi_tax_id: int = 9606,
        score_threshold: Union[
            Number,
            Literal[
                'highest_confidence',
                'high_confidence',
                'medium_confidence',
                'low_confidence',
            ]
        ] = 'highest_confidence',
        physical_interaction_score: bool = True,
    ) -> List[tuple]:
    """
    Downloads chemical-protein links (with detailed subscores).
    The combined physical interaction score is defined
    between the interactors for which we have evidence of their binding.

    Args
        score_threshold: Minimum required interaction score. user can use
            pre-defined confidence limits or can define a custom value.
    """

    confidence = {
        'highest_confidence': 900,
        'high_confidence': 700,
        'medium_confidence': 400,
        'low_confidence': .150,
    }

    min_score = confidence.get(score_threshold, score_threshold)

    StitchLinksInteraction = collections.namedtuple(
        'StitchLinksInteraction',
        (
            'partner_a',
            'partner_b',
            'experimental',
            'prediction',
            'database',
            'textmining',
            'combined_score',
            'physical_combined_score',
        ),
    )

    if physical_interaction_score:

        phy_links = dict(
            (
                (s.partner_a, s.partner_b),
                s.score
            )
            for s in stitch_actions_interactions()
            if s.mechanism == 'binding'
        )

    url = urls.urls['stitch']['links'] % ncbi_tax_id
    c = curl.Curl(url, silent = False, large = True)
    _ = next(c.result)

    sep = re.compile(r'[sm\.]')

    for l in c.result:

        if hasattr(l, 'decode'):

            l = l.decode('utf-8')

        l = l.strip().split('\t')

        if int(l[6]) < min_score:

            continue

        a = sep.split(l[0])[1]
        b = sep.split(l[1])[1]

        phy_score = (
            phy_links.get((a,b), phy_links.get((b, a), None))
                if physical_interaction_score else
            None
        )

        yield StitchLinksInteraction(
            partner_a = a,
            partner_b = b,
            experimental = int(l[2]),
            prediction = int(l[3]),
            database = int(l[4]),
            textmining = int(l[5]),
            combined_score = int(l[6]),
            physical_combined_score = phy_score,
        )
