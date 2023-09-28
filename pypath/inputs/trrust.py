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

import collections

from pypath.share import curl
from pypath.resources.urls import urls
from pypath.share import session
from pypath.utils import taxonomy

_log = session.Logger(name = 'trrust_input')._log


def trrust_interactions(
        organism: str | int = 'human',
    ) -> list[TrrustInteraction]:
    """
    Gene regulatory interactions from the TRRUST v2 database.

    https://academic.oup.com/nar/article/46/D1/D380/4566018

    Args:
        organism:
            Name or NCBI Taxonomy ID of the organism. Human and mouse
            are available in TRRUST.
    """

    organisms = {'human', 'mouse'}
    _organism = taxonomy.ensure_common_name(organism, lower = True)

    if _organism not in organisms:

        err = f'Only human and mouse are availble in TRRUST, not `{organism}`.'
        _log(err)
        raise ValueError(err)

    class TrrustInteraction(
            collections.namedtuple(
            'TrrustInteractionBase',
            ('source_genesymbol', 'target_genesymbol', 'effect', 'references'),
        )
    ):

        def __new__(cls, line):

            line = line.strip('\n ').split('\t')
            refs = tuple(sorted(line[-1].split(';')))

            return super().__new__(cls, *line[:-1], refs)


    url = urls['trrust']['tsv_url'] % _organism

    c = curl.Curl(
        url,
        silent = False,
        large = True,
        encoding = 'utf-8',
        default_mode = 'r',
    )

    interactions = [TrrustInteraction(line) for line in c.result]

    return interactions


def trrust_human():

    return trrust_interactions('human')


def trrust_mouse():

    return trrust_interactions('mouse')
