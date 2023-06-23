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

"""
CollecTRI is a comprehensive resource of TF-target interactions.
"""

from typing import Literal

import re
import collections
import itertools

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.internals.intera as intera
import pypath.utils.mapping as mapping
import pypath.share.session as session

_log = session.Logger(name = 'collectri_input')._log

# Based on literature by Sophia
# https://www.sciencedirect.com/science/article/abs/pii/S0304419X19300526
# https://www.nature.com/articles/1209933
COMPLEXES = {
    'AP1': {
        'JUN-FOS',
        'JUNB-FOS',
        'JUND-FOS',
        'JUN-FOSB',
        'JUNB-FOSB',
        'JUND-FOSB',
        'JUN-FOSL1',
        'JUNB-FOSL1',
        'JUND-FOSL1',
        'JUN-FOSL2',
        'JUNB-FOSL2',
        'JUND-FOSL2',
        'JUN-JUN',
        'JUN-JUNB',
        'JUN-JUND',
        'JUNB-JUNB',
        'JUNB-JUND',
        'JUND-JUND',
    },
    'NFKB': {
        'RELA-RELA',
        'RELA-REL',
        'RELA-RELB',
        'RELA-NFKB1',
        'RELA-NFKB2',
        'REL-REL',
        'REL-RELB',
        'REL-NFKB1',
        'REL-NFKB2',
        'RELB-RELB',
        'RELB-NFKB1',
        'RELB-NFKB2',
        'NFKB1-NFKB1',
        'NFKB1-NFKB2',
        'NFKB2-NFKB2',
    },
}


def collectri_raw(
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

    CollectriRecord = collections.namedtuple(
        'CollectriRecord',
        (
            'tf',
            'target',
            'effect',
            'tf_category',
            'resources',
            'pubmed',
            'sign_decision',
            'target_type',
        ),
    )

    url = urls.urls['collectri']['url']
    c = curl.Curl(
        url,
        silent = False,
        large = True,
    )
    _ = next(c.result)

    result = []

    for l in c.result:

        l = l.strip().split(',')
        mmirna = remirna.match(l[1])

        if (
            (mmirna and not mirna) or
            (not mmirna and not protein_coding)
        ):

            continue

        target_id = f'hsa-miR-{mmirna.group(1)}' if mmirna else l[1]

        result.append(
            CollectriRecord(
                tf = l[0],
                target = target_id,
                effect = int(l[2]),
                tf_category = l[3],
                resources = l[4].replace('DoRothEA_A', 'DoRothEA-A'),
                pubmed = l[5],
                sign_decision = l[6],
                target_type = 'mirna' if mmirna else 'protein',
            )
        )

    return result


def collectri_interactions(
        protein_coding: bool = True,
        mirna: bool = False,
    ) -> list[tuple]:
    """
    TF-target interactions from the CollecTRI database.

    While `collectri_raw` returns the records in the same format as in the
    original data, here we translate identifiers to UniProt IDs, and use
    `Complex` objects to represent protein complexes.

    Args:
        protein_coding:
            Include regulation of protein coding genes.
        mirna:
            Include regulation of miRNA coding genes.
    """

    CollectriInteraction = collections.namedtuple(
        'CollectriInteraction',
        (
            'tf',
            'target',
            'effect',
            'tf_category',
            'resources',
            'pubmed',
            'sign_decision',
            'target_type',
        ),
    )


    def process_complex(name):

        result = []

        for var in COMPLEXES[name]:

            uniprots = [
                mapping.map_name(comp, 'genesymbol', 'uniprot')
                for comp in var.split('-')
            ]

            if all(uniprots):

                result.extend(list(itertools.product(*uniprots)))

            else:

                _log(
                    'Failed to translate all components of '
                    f'complex `{name}` (components: {var}).'
                )

        return set(result)


    for rec in collectri_raw(protein_coding = protein_coding, mirna = mirna):

        tf_uniprots = (
            process_complex(rec.tf)
                if rec.tf in COMPLEXES else
            mapping.map_name(rec.tf, 'genesymbol', 'uniprot')
        )

        target_uniprots = (
            (rec.target,)
                if rec.target_type == 'mirna' else
            mapping.map_name(rec.target, 'genesymbol', 'uniprot')
        )

        for tf_u, t_u in itertools.product(tf_uniprots, target_uniprots):

            if isinstance(tf_u, tuple):

                tf_u = intera.Complex(
                    components = tf_u,
                    sources = 'CollecTRI',
                )

            yield CollectriInteraction(tf_u, t_u, *rec[2:])

