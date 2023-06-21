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

from __future__ import annotations

from typing import Iterable, Literal

import collections

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.common as common


CONFIDENCE_THRESHOLDS = {
    'highest_confidence': 900,
    'high_confidence': 700,
    'medium_confidence': 400,
    'low_confidence': .150,
}


def string_effects(
        ncbi_tax_id: int = 9606,
        stimulation: str | Iterable[str] = 'activation',
        inhibition: str | Iterable[str] = 'inhibition',
        exclude: str | Iterable[str] = 'expression',
        score_threshold: int = 0,
    ) -> list[tuple]:

    StringEffectsInteraction = collections.namedtuple(
        'StringEffectsInteraction',
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
                    StringEffectsInteraction(
                        source = source,
                        target = target,
                        effect = effect,
                    )
                )

    return effects


def string_links_interactions(
        ncbi_tax_id: int = 9606,
        score_threshold: int | Literal[
            'highest_confidence',
            'high_confidence',
            'medium_confidence',
            'low_confidence',
            ] = 'highest_confidence',
        physical_interaction_score: bool = True,
    ) -> list[tuple]:
    """
    Downloads protein network data, including subscores per channel.
    The output contains both functional and physical protein associations.
    The combined physical interaction score is defined between the proteins
    for which we have evidence of their binding or forming a physical complex.

    Args
        score_threshold: Minimum required interaction score. user can use
            pre-defined confidence limits or can define a custom value.
    """

    StringLinksInteraction = collections.namedtuple(
        'StringLinksInteraction',
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

        phy_links = dict(
            (
                (i.protein_a,i.protein_b),
                i.combined_score
            )
            for i in
            string_physical_interactions(
                ncbi_tax_id = ncbi_tax_id,
                score_threshold = 0,
            )
        )

    url = urls.urls['string']['links'] % ncbi_tax_id
    c = curl.Curl(url, silent = False, large = True)
    _ = next(c.result)

    min_score = CONFIDENCE_THRESHOLDS.get(score_threshold, score_threshold)

    for l in c.result:

        l = l.strip().split(' ')
        prot_a_id = l[0].split('.')[1]
        prot_b_id = l[1].split('.')[1]

        if int(l[9]) < min_score:

            continue

        phy_score = (
            phy_links.get((prot_a_id, prot_b_id), None)
                if physical_interaction_score else
            None
        )

        yield StringLinksInteraction(
            protein_a = prot_a_id,
            protein_b = prot_b_id,
            neighborhood_score = int(l[2]),
            fusion = int(l[3]),
            cooccurence = int(l[4]),
            coexpression = int(l[5]),
            experimental = int(l[6]),
            database = int(l[7]),
            textmining = int(l[8]),
            combined_score = int(l[9]),
            physical_combined_score = phy_score,
        )


def string_physical_interactions(
        ncbi_tax_id: int = 9606,
        score_threshold: int | Literal[
            'highest_confidence',
            'high_confidence',
            'medium_confidence',
            'low_confidence',
        ] = 'highest_confidence',
    ) -> list[tuple]:
    """
    Downloads protein physical subnetwork data, including subscores per
    channel. The interactions indicate that the proteins are part of a
    physical complex.

    Args
        score_threshold: Minimum required interaction score. user can use
            pre-defined confidence limits or can define a custom value.
    """

    StringPhysicalInteraction = collections.namedtuple(
        'StringPhysicalInteraction',
        (
            'protein_a',
            'protein_b',
            'experimental',
            'database',
            'textmining',
            'combined_score',
        ),
    )

    links = []

    url = urls.urls['string']['physical_links'] % ncbi_tax_id
    c = curl.Curl(url, silent = False, large = True)
    _ = next(c.result)

    min_score = CONFIDENCE_THRESHOLDS.get(score_threshold, score_threshold)

    for l in c.result:

        l = l.strip().split(' ')

        if int(l[5]) >= min_score:

            links.append(
                StringPhysicalInteraction(
                    protein_a= l[0].split('.')[1],
                    protein_b= l[1].split('.')[1],
                    experimental= int(l[2]),
                    database= int(l[3]),
                    textmining= int(l[4]),
                    combined_score= int(l[5]),
                )
            )

    return links


def string_species() -> dict[int, str]:
    """
    Downloads list of organisms in STRING.

    Returns
        Dict of tax ids as keys and scientific names of organisms as values.
    """

    species = {}

    url = urls.urls['string']['species']
    c = curl.Curl(url, silent = False, large = True)
    _ = next(c.result)

    for l in c.result:

        l = l.strip().split('\t')
        tax_id = l[0]
        official_name = l[3]

        species[tax_id] = official_name

    return species