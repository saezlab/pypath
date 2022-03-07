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

import re
import collections

import pypath.resources.urls as urls
import pypath.share.curl as curl


def stitch_interactions(threshold = None):

    StitchInteraction = collections.namedtuple(
        'StitchInteraction',
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

        try:
            a = sep.split(l[0])[1]
            b = sep.split(l[1])[1]

        except IndexError:
            print(l[1])

        if l[4] == 'f':
            a, b = b, a

        yield StitchInteraction(
            partner_a = a,
            partner_b = b,
            mechanism = l[2],
            action = l[3],
            score = int(l[5]),
        )

def stitch_links(
    ncbi_tax_id = 9606,
    score_threshold = 'highest_confidence',
    physical_interaction_score= True,
):

    """
    Downloads chemical-protein links (with detailed subscores).  
    The combined physical interaction score is defined 
    between the interactors for which we have evidence of their binding.

    :param int,str score_threshold: 
        minimum required interaction score. user can use pre-defined confidence limits or can define a custom value.


    """

    confidence= {'highest_confidence':900,
        'high_confidence':700,
        'medium_confidence':400,
        'low_confidence':0.150}

    if score_threshold not in confidence:
        confidence[score_threshold]= score_threshold

    StitchLinks = collections.namedtuple(
        'StitchLinks',
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
        phy_links=[s for s in stitch_interactions() if s.mechanism =="binding"]
        phy_links_dict={}
        for i in phy_links:
            phy_links_dict[(i.partner_a, i.partner_b)] = i.score

    url = urls.urls['stitch']['links'] % ncbi_tax_id
    c = curl.Curl(url, silent = False, large = True)
    _ = next(c.result)

    sep = re.compile(r'[sm\.]')
    

    for l in c.result:
        if hasattr(l, 'decode'):
            l = l.decode('utf-8')
        l = l.strip().split('\t')

        if int(l[6]) < confidence[score_threshold]:
            continue
        try:
            a = sep.split(l[0])[1]
            b = sep.split(l[1])[1]

        except IndexError:
            print(l[1])

        if not physical_interaction_score:
            phy_score='' 
        else:
            a,b = b,a
            if (a,b) in phy_links_dict:
                phy_score= phy_links_dict[(a,b)]
            else:
                phy_score=''

        yield StitchLinks(
            partner_a= a,
            partner_b= b,
            experimental= int(l[2]),
            prediction= int(l[3]),
            database= int(l[4]),
            textmining= int(l[5]),
            combined_score= int(l[6]),
            physical_combined_score= phy_score,
                )
