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
Structures from the Human Metabolome Database (HMDB).
"""

import pypath.formats.sdf as sdfparser
import pypath.resources.urls as urls
import pypath.share.curl as curl


def sdf():
    """
    Download and open the SDF file with all HMDB structures.
    """

    url = urls.urls['hmdb']['sdf']
    c = curl.Curl(
        url,
        large = True,
        silent = False,
        default_mode = 'rb',
        files_needed = ['structures.sdf'],
    )

    return sdfparser.SdfReader(
        c.result['structures.sdf'],
        names = {
            'DATABASE_ID': 'id',
            'SMILES': 'smiles',
            'GENERIC_NAME': 'generic_name',
        }
    )
