#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2023
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: Dénes Türei (turei.denes@gmail.com)
#           Nicolàs Palacio
#           Sebastian Lobentanzer
#           Erva Ulusoy
#           Olga Ivanova
#           Ahmet Rifaioglu
#           Sophia Müller-Dott
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

from __future__ import annotations

"""
CollecTRI is a comprehensive resource of TF-target interactions.
"""

from typing import Literal

import re
import collections

import pypath.resources.urls as urls
import pypath.share.curl as curl


def collectri_interactions(
        protein_coding: bool = True,
        mirna: bool = False,
    ) -> list[tuple]:
    """
    TF-target interactions from the CollecTRI database.

    Args:
        protein_coding:
            Include regulation of protein coding genes.
        mirna:
            Include regulation of miRNA coding genes.
    """


    remirna = re.compile('^MIR(\d+)$')

    CollectriInteraction = collections.namedtuple(
        'CollectriInteraction',
        (
            'tf',
            'target',
            'effect',
            'TF_category',
            'resources',
            'pubmed',
            'target_type',
        ),
    )

    url = urls.urls['collectri']['url']
    c = curl.Curl(
        url,
        silent = False,
        large = True,
    )

    result = []

    for l in c.result:

        l = l.strip().split(',')
        mmirna = remirna.match(l[1])

        if (
            (mmirna and not mirna) or
            (not mmirna and not protein)
        ):

            continue

        target_id = f'hsa-miR-{mmirna.group(1)}' if mmirna else l[1]

        result.append(
            CollectriInteraction(
                tf = l[0],
                target = target_id,
                effect = int(l[2]),
                tf_category = l[3],
                resources = l[4],
                pubmed = l[5],
                target_type = 'mirna' if mmirna else 'protein',
            )
        )

    return result
