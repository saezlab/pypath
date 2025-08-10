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

from pypath.resources import urls
from pypath.share import curl
from pypath.utils import taxonomy

HintInteraction = collections.namedtuple(
    'HintInteraction',
    [
        'id_a',
        'id_b',
        'pmids',
        'quality',
    ],
)


def hint_raw(organism = 9606):

    organism_latin = taxonomy.ensure_latin_name(organism)

    if organism_latin is None:

        raise ValueError('Unknown organism: %s' % organism)

    organism_hint = organism_latin.title().replace(' ', '')
    url = urls.urls['hint']['url'] % organism_hint
    c = curl.Curl(url = url, large = True, silent = False)
    _ = next(c.result)

    for line in c.result:

        yield line.strip().split('\t')


def _hint_parse_refs(quality, refs):

    return [
        this_ref[0]
        for ref in refs
        if (this_ref := ref.split(':'))[2] == quality
    ]


def hint_interactions(organism = 9606):

    for raw in hint_raw(organism = organism):

        refs = raw[4].split('|')

        for quality in ('LC', 'HT'):

            this_refs = _hint_parse_refs('LC', refs)

            if not this_refs:

                continue

            yield HintInteraction(
                id_a = raw[0],
                id_b = raw[1],
                pmids = this_refs,
                quality = quality,
            )
