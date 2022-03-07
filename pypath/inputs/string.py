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
import pypath.share.common as common


def string_effects(
    ncbi_tax_id = 9606,
    stimulation = ('activation',),
    inhibition = ('inhibition',),
    exclude = ('expression',),
    score_threshold = 0,
):

    StringInteraction = collections.namedtuple(
        'StringInteraction',
        (
            'source',
            'target',
            'effect',
        ),
    )

    effects = []

    stimulation = common.to_set(stimulation)
    inhibition = common.to_set(inhibition)
    exclude = common.to_set(exclude)

    url = urls.urls['string']['actions'] % ncbi_tax_id
    c = curl.Curl(url, silent = False, large = True)
    _ = next(c.result)

    for l in c.result:

        if hasattr(l, 'decode'):

            l = l.decode('ascii')

        l = l.strip().split('\t')

        if l and l[4] == 't' and int(l[6]) >= score_threshold:

            effect = (
                '+'
                    if l[2] in stimulation else
                '-'
                    if l[2] in inhibition else
                '*'
                    if l[2] not in exclude else
                None
            )

            source = l[0].split('.')[1] if l[5] == 't' else l[1].split('.')[1]
            target = l[1].split('.')[1] if l[5] == 't' else l[0].split('.')[1]

            if effect is not None:

                effects.append(
                    StringInteraction(
                        source = source,
                        target = target,
                        effect = effect,
                    )
                )

    return effects


def string_links(
    ncbi_tax_id = 9606,
    score_threshold = 'highest_confidence',
    physical_interaction_score= True,
):

    """
    Downloads protein network data, including subscores per channel.  
    The output contains both functional and physical protein associations.
    The combined physical interaction score is defined 
    between the proteins for which we have evidence of their binding or forming a physical complex.
    
    :param int,str score_threshold: 
        minimum required interaction score. user can use pre-defined confidence limits or can define a custom value.
    """

    confidence= {'highest_confidence':900,
        'high_confidence':700,
        'medium_confidence':400,
        'low_confidence':0.150}

    StringLinks = collections.namedtuple(
        'StringLinks',
        (
            'protein_a',
            'protein_b',
            'neighborhood_score',
            'fusion',
            'cooccurence',
            'coexpression',
            'experimental',
            'database',
            'textmining',
            'combined_score',
            'physical_combined_score',
        ),
    )
    if physical_interaction_score:
        phy_links=string_physical_links(ncbi_tax_id, score_threshold=0)
        phy_links_dict={}
        for i in phy_links:
            phy_links_dict[(i.protein_a,i.protein_b)] = i.combined_score

    url = urls.urls['string']['links'] % ncbi_tax_id
    c = curl.Curl(url, silent = False, large = True)
    _ = next(c.result)

    
    if score_threshold not in confidence:
        confidence[score_threshold]= score_threshold

    for l in c.result:
        l = l.strip().split(' ')
        prot_a_id=l[0].split('.')[1]
        prot_b_id=l[1].split('.')[1]
        if int(l[9]) < confidence[score_threshold]:
            continue

        if physical_interaction_score and (prot_a_id,prot_b_id) in phy_links_dict:
            phy_score = phy_links_dict[(prot_a_id,prot_b_id)]
        else:
            phy_score=''


        yield StringLinks(
            protein_a= prot_a_id,
            protein_b= prot_b_id,
            neighborhood_score= int(l[2]),
            fusion= int(l[3]),
            cooccurence= int(l[4]),
            coexpression= int(l[5]),
            experimental= int(l[6]),
            database= int(l[7]),
            textmining= int(l[8]),
            combined_score= int(l[9]),
            physical_combined_score=phy_score,
        )


def string_physical_links(
    ncbi_tax_id = 9606,
    score_threshold = 'highest_confidence'
):

    """
    Downloads protein physical subnetwork data, including subscores per channel.
    The interactions indicate that the proteins are part of a physical complex.
    
    :param int,str score_threshold: 
        minimum required interaction score. user can use pre-defined confidence limits or can define a custom value.
    """

    confidence= {'highest_confidence':900,
        'high_confidence':700,
        'medium_confidence':400,
        'low_confidence':0.150}

    StringLinks = collections.namedtuple(
        'StringLinks',
        (
            'protein_a',
            'protein_b',
            'experimental',
            'database',
            'textmining',
            'combined_score'
        ),
    )

    links=[]


    url = urls.urls['string']['physical_links'] % ncbi_tax_id
    c = curl.Curl(url, silent = False, large = True)
    _ = next(c.result)
    
    if score_threshold in confidence:
        for l in c.result:
            l = l.strip().split(' ')
            if int(l[5]) >= confidence[score_threshold]:
                links.append(
                    StringLinks(
                        protein_a= l[0].split('.')[1],
                        protein_b= l[1].split('.')[1],
                        experimental= int(l[2]),
                        database= int(l[3]),
                        textmining= int(l[4]),
                        combined_score= int(l[5])
                        )
                    )

    else:
        confidence["custom"]= score_threshold
        for l in c.result:
            l = l.strip().split(' ')
            if int(l[5]) >= confidence["custom"]:
                links.append(
                    StringLinks(
                        protein_a= l[0].split('.')[1],
                        protein_b= l[1].split('.')[1],
                        experimental= int(l[2]),
                        database= int(l[3]),
                        textmining= int(l[4]),
                        combined_score= int(l[5])
                        )
                    )
    return links


