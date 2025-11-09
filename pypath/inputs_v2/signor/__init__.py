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

"""
Download and process data from SIGNOR database.

This v2 module emits Entity records according to the schema defined in
pypath.internals.silver_schema.
"""

from ._parsers import (
    signor_complexes,
    signor_protein_families,
    signor_phenotypes,
    signor_stimuli,
    signor_interactions
)

__all__ = [
    'SIGNOR_RESOURCE',
    'signor_interactions',
    'signor_complexes',
    'signor_protein_families',
    'signor_phenotypes',
    'signor_stimuli',
]
