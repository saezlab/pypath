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

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.utils.mapping as mapping


def mcam_cell_adhesion_molecules():
    
    url = urls.urls['mcam']['url']
    
    c = curl.Curl(url, silent = False, large = True)
    
    _ = next(c.result)
    
    return {
        uniprot
        for line in c.result
        for _uniprot in line.split('\t', maxsplit = 1)[0].split(';')
        for uniprot in mapping.map_name(
            _uniprot.strip(),
            'uniprot',
            'uniprot',
        )
    }
