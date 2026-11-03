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

import zipfile

import pypath.formats.sdf as sdfparser
import pypath.resources.urls as urls
from pypath.share.downloads import dm


__all__ = [
    'lipidmaps_sdf',
    'lipidmaps_raw',
]


def lipidmaps_sdf():
    """
    Download and open the SDF file with all LipidMaps (LMSD) structures.
    """

    url = urls.urls['lipidmaps']['lmsd']

    # Download file using download manager
    file_path = dm.download(
        url,
        filename='structures.zip',
        subfolder='lipidmaps',
    )

    # Extract and read the SDF file from the zip
    with zipfile.ZipFile(file_path) as zf:
        sdf_file = zf.open('structures.sdf')
        return sdfparser.SdfReader(
            sdf_file,
            names = {
                'HMDB_ID': 'hmdb_id',
                'PUBCHEM_CID': 'pubchem',
                'SWISSLIPIDS_ID': 'swisslipids',
                'LM_ID': 'lipidmaps',
                'ABBREVIATION': 'abbreviation',
            }
        )


def lipidmaps_raw():

    yield from lipidmaps_sdf()
