#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright 2014-2024
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

"""
Named tuple definitions for Metabolic Atlas data structures.
"""

from __future__ import annotations

import collections

__all__ = [
    'MetatlasModel',
    'MetatlasGem',
    'GemReaction',
    'GemMetabolite',
]


MetatlasModel = collections.namedtuple(
    'MetatlasModel',
    [
        'id',
        'name',
        'short_name',
        'organism',
        'tissue',
        'cell_type',
        'condition',
        'year',
        'reaction_count',
        'metabolite_count',
        'gene_count',
        'files',
    ],
)


MetatlasGem = collections.namedtuple(
    'MetatlasGem',
    [
        'model_name',
        'git_host',
        'git_repo',
        'git_url',
        'latest_version',
        'commits',
        'contributors',
    ],
)


GemReaction = collections.namedtuple(
    'GemReaction',
    [
        'id',
        'name',
        'metabolites',
        'lower_bound',
        'upper_bound',
        'gene_reaction_rule',
        'subsystem',
        'eccodes',
    ],
)


GemMetabolite = collections.namedtuple(
    'GemMetabolite',
    [
        'id',
        'name',
        'compartment',
        'formula',
        'charge',
    ],
)


