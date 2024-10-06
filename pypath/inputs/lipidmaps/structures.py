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
Structures from the LIPID MAPS Structure Database (LMSD).
"""

import pypath.formats.sdf as sdfparser
import pypath.resources.urls as urls
import pypath.share.curl as curl


def sdf():
    """
    Download and open the SDF file with all LipidMaps (LMSD) structures.
    """

    url = urls.urls['lipidmaps']['lmsd']
    c = curl.Curl(
        url,
        large = True,
        silent = False,
        default_mode = 'rb',
        compr = 'zip',
        files_needed = ['structures.sdf'],
    )

    return sdfparser.SdfReader(
        c.result['structures.sdf'],
        names = {
            'HMDB_ID': 'hmdb_id',
            'PUBCHEM_CID': 'pubchem',
            'SWISSLIPIDS_ID': 'swisslipids',
            'LM_ID': 'lipidmaps',
            'ABBREVIATION': 'abbreviation',
        }
    )
