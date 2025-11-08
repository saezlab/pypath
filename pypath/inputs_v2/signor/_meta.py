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
Metadata for SIGNOR database resource.

This module defines the Resource metadata for SIGNOR according to the
schema defined in pypath.internals.silver_schema.
"""

from pypath.internals.silver_schema import Resource
from omnipath_build.utils.cv_terms import LicenseCV, UpdateCategoryCV

__all__ = ['SIGNOR_RESOURCE']


SIGNOR_RESOURCE = Resource(
    id='signor',
    name='SIGNOR',
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    publication='PMID:31665520',  # SIGNOR 2.0 paper
    url='https://signor.uniroma2.it/',
    description=(
        'SIGNOR (SIGnaling Network Open Resource) is a comprehensive '
        'resource of causal relationships between biological entities '
        'with a focus on signaling pathways. It provides manually curated '
        'interactions with mechanistic details including protein-protein '
        'interactions, post-translational modifications, transcriptional '
        'regulation, and small molecule effects.'
    ),
)
