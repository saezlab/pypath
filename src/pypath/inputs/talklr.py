#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2021
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

import pyreadr


import pypath.share.curl as curl
import pypath.resources.urls as urls


def talklr_raw():

    url = urls.urls['talklr']['url']
    c = curl.Curl(url, large = True, silent = False)
    rdata_path = c.fileobj.name
    c.fileobj.close()

    rdata = pyreadr.read_r(rdata_path)['receptor_ligand']
    rdata.columns = [col.replace('.', '_') for col in rdata.columns]

    return rdata

