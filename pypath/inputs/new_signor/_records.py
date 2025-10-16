#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright 2014-2025
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

__all__ = [
    'SignorComplex',
    'SignorProteinFamily',
    'SignorPhenotype',
    'SignorStimulus',
]


SignorComplex = collections.namedtuple(
    'SignorComplex',
    ['complex_id', 'name', 'components']
)

SignorProteinFamily = collections.namedtuple(
    'SignorProteinFamily',
    ['family_id', 'name', 'members']
)

SignorPhenotype = collections.namedtuple(
    'SignorPhenotype',
    ['phenotype_id', 'name', 'description']
)

SignorStimulus = collections.namedtuple(
    'SignorStimulus',
    ['stimulus_id', 'name', 'description']
)
