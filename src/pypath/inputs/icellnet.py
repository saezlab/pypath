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

import collections

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.inputs.common as inputs_common


IcellnetRecord = collections.namedtuple(
    'IcellnetRecord',
    [
        
    ]
)


def icellnet():
    
    url = urls.urls['icellnet']['url']
    
    c = curl.Curl(url, silent = False, large = True)
    
    xls = c.fileobj
    xlsfile = xls.name
    xls.close()
    tbl = inputs_common.read_xls(xlsfile)
    
    return filter(lambda l: len(l[-1]) > 0, map(lambda l: l[:7], tbl[2:]))
