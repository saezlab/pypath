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

__all__ = [
    'NAME_AUTHORITIES',
    'POSITIVE_REGULATION',
    'NEGATIVE_REGULATION',
]


NAME_AUTHORITIES = {
    'Human': 'HGNC',
    'Rat': 'RGD',
    'Mouse': 'MGI',
}

POSITIVE_REGULATION = {
    'agonist',
    'activator',
    'partial agonist',
    'inverse antagonist',
    'full agonist',
    'activation',
    'irreversible agonist',
    'positive',
    'biased agonist',
}

NEGATIVE_REGULATION = {
    'inhibitor',
    'antagonist',
    'competitive',
    'feedback inhibition',
    'inhibition',
    'irreversible inhibition',
    'inverse agonist',
    'negative',
    'weak inhibition',
    'reversible inhibition',
    'voltage-dependent inhibition',
    'pore blocker',
}
