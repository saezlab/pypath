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

import collections

__all__ = [
    'LipidmapsLipid',
]


LipidmapsLipid = collections.namedtuple(
    'LipidmapsLipid',
    [
        'id',
        'name',
        'abbreviation',
        'formula',
        'exact_mass',
        'category',
        'main_class',
        'smiles',
        'inchi',
        'inchikey',
        'synonyms',
        'iupac',
        'pubchem',
        'chebi',
    ]
)
