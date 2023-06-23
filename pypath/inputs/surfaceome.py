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

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.utils.mapping as mapping
import pypath.inputs.common as inputs_common


def surfaceome_annotations():
    """
    Downloads the "In silico human surfaceome".
    Dict with UniProt IDs as key and tuples of surface prediction score,
    class and subclass as values (columns B, N, S and T of table S3).
    """

    url = urls.urls['surfaceome']['url']
    c = curl.Curl(url, large = True, silent = False)
    xlsname = c.fname
    del(c)
    raw = inputs_common.read_xls(xlsname, 'in silico surfaceome only')[2:]

    return dict(
        (
            uniprot, # uniprot
            (
                float(r[13]), # score
                r[18] if r[18] else None, # class
                set(r[19].replace('KInase', 'Kinase').split(';'))
                    if r[19] else
                set(), # subclass
            )
        )
        for r in raw
        for uniprot in mapping.map_name(r[1], 'uniprot', 'uniprot')
    )
