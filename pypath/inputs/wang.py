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

import collections

from typing import List, Literal, Union

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.settings as settings
import pypath.inputs.embopress as embo
import pypath.inputs.ca1 as ca1
import pypath.utils.mapping as mapping


KEY = {
    'NA': None,
    'Neucleus': 'Nucleus',
    'Ribosomes': 'Ribosome',
    'Vesicles': 'Vesicle',
    'Endoplasmic reticulum': 'ER',
    'Not available': None,
    'Not Available': None,
    'Mitochondrial': 'Mitochondrion',
    'Mitochondria': 'Mitochondrion',
    'Anti-Apoptic': 'Anti-apoptotic',
}


ca1_interactions = ca1.ca1_interactions


def hsn_interactions(
        source: Literal['rescued', 'researchgate'] = 'rescued',
    ) -> List[tuple]:
    """
    Downloads and processes HumanSignalingNetwork version 6
    (published 2014 Jan by Edwin Wang).

    Args
        source: The same file is available from two domains: the OmniPath
            rescued repository and Research Gate. These both are secondary
            sources, the dataset is not available any more from its original
            site of publication, which was the old webpage of the Wang Lab.

    Details:
        This dataset is identical to the one returned by `wang_interactions`,
        but it does not contain function and localization details.
    """

    effects = {
        'Pos': '+',
        'Neg': '-',
        'Phy': '0',
    }

    class HsnInteraction(
            collections.namedtuple(
                'HsnInteractionBase',
                (
                    'genesymbol_source',
                    'genesymbol_target',
                    'entrez_source',
                    'entrez_target',
                    'effect',
                ),
            ),
        ):

        def __new__(cls, *args):

            args = args[0] if len(args) == 1 else args
            identifiers = dict(zip(
                (
                    'entrez_source',
                    'genesymbol_source',
                    'entrez_target',
                    'genesymbol_target',
                ),
                args[:-1]
            ))

            return super(HsnInteraction, cls).__new__(
                cls,
                **identifiers,
                effect = effects.get(args[-1], args[-1]),
            )


    url = urls.urls['hsn'][source]
    c = curl.Curl(
        url,
        silent = False,
        large = True,
        req_headers = [settings.get('user_agent')],
    )
    _ = next(c.result)
    result = [HsnInteraction(*r.strip().split(',')) for r in c.result if r]

    return result


def wang_interactions() -> List[tuple]:
    """
    Downloads and processes Wang Lab HumanSignalingNetwork.
    Returns list of interactions as tuples of source, target and effect.
    """

    url = urls.urls['wang']['rescued']
    c = curl.Curl(url, silent = False)
    data = c.result
    data = data.split('\n')

    return _wang_process(data)


def cui_interactions() -> List[tuple]:
    """
    Interactions from Supplementary Table 9 of Cui et al. 2007
    (https://www.embopress.org/doi/full/10.1038/msb4100200).
    """

    raw = embo.embopress_supplementary(
        url = urls.urls['wang']['cui'],
        init_url = urls.urls['wang']['cui_init'],
        sheet = 'Supplementary Table 9',
    )

    return _wang_process(raw)


def _wang_process(raw: List[Union[List, str]]) -> List[tuple]:

    Node = collections.namedtuple(
        'Node',
        (
            'genesymbol',
            'entrez',
            'function',
            'location',
        ),
    )

    WangInteraction = collections.namedtuple(
        'WangInteraction',
        (
            'genesymbol_source',
            'genesymbol_target',
            'entrez_source',
            'entrez_target',
            'effect',
            'function_source',
            'location_source',
            'function_target',
            'location_target',
        ),
    )

    key = KEY.copy()

    result = []
    nodes = {}
    reading_nodes = False
    reading_edges = False
    reading_key = False

    _key = lambda y: (lambda x: key.get(x, x))(key.get(y, y))

    for l in raw:

        if not l or not l[0] or (hasattr(l, 'strip') and not l.strip()):

            reading_key = False
            reading_nodes = False
            reading_edges = False

        l = l.split(',') if hasattr(l, 'split') else l

        if reading_key:

            key[l[0]] = l[1]

        elif reading_nodes:

            nodes[l[0]] = Node(
                genesymbol = l[1],
                entrez = l[2].split('.')[0],
                function = _key(l[3]),
                location = _key(l[4]),
            )

        elif reading_edges:

            src = nodes[l[0]]
            tgt = nodes[l[1]]

            result.append(
                WangInteraction(
                    genesymbol_source = src.genesymbol,
                    genesymbol_target = tgt.genesymbol,
                    entrez_source = src.entrez,
                    entrez_target = tgt.entrez,
                    effect = l[2].replace('_', '-').split('.')[0],
                    function_source = src.function,
                    location_source = src.location,
                    function_target = tgt.function,
                    location_target = tgt.location,
                )
            )

        if l[0].startswith('Node'):

            reading_key = False
            reading_nodes = True

        if l[0].startswith('From'):

            reading_key = False
            reading_nodes = False
            reading_edges = True

        if l[0].startswith('Notes:'):

            reading_key = True

    return result


def wang_annotations():

    NA_VALUES = {None, 'NA', 'Not available'}

    record = collections.namedtuple(
        'WangAnnotation',
        (
            'function',
            'location',
        ),
    )

    result = collections.defaultdict(set)

    for dataset in ('ca1', 'cui', 'wang'):

        func = globals()['%s_interactions' % dataset]
        data = func()

        for i in data:

            for side in ('source', 'target'):

                a_uniprot = 'uniprot_%s' % side
                a_entrez = 'entrez_%s' % side
                a_function = 'function_%s' % side
                a_location = 'location_%s' % side

                function = getattr(i, a_function)
                location = getattr(i, a_location)

                if function in NA_VALUES and location in NA_VALUES:

                    continue

                uniprots = (
                    (getattr(i, a_uniprot),)
                        if hasattr(i, a_uniprot) else
                    mapping.map_name(getattr(i, a_entrez), 'entrez', 'uniprot')
                )

                for uniprot in uniprots:

                    result[uniprot].add(
                        record(
                            function = function,
                            location = location,
                        )
                    )

    return dict(result)
