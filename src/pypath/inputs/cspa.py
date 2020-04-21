#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2020
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#                  Nicolàs Palacio
#                  Olga Ivanova
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.inputs.common as inputs_common
import pypath.utils.mapping as mapping
import pypath.utils.taxonomy as taxonomy


def get_cspa(organism = 9606):

    sheets = {
        'Human': 'Table A',
        'Mouse': 'Table B',
    }

    str_organism = taxonomy.taxids[organism].capitalize()

    url = urls.urls['cspa']['url']
    c = curl.Curl(url, large = True, silent = False)
    xlsname = c.fname
    del(c)
    raw = inputs_common.read_xls(xlsname, sheets[str_organism])[1:]

    return mapping.map_names((r[1] for r in raw), 'uniprot', 'uniprot')
