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
import pypath.inputs.common as inputs_common


def _slc_raw(table = 1, sheet = 1):

    url = urls.urls['slc'][f'table_s{table}']
    c = curl.Curl(url, large = True, silent = False)
    content = inputs_common.read_xls(c.fileobj.name, sheet = sheet)

    return content


def slc_annotation():

    return _slc_raw(table = 1, sheet = 1)


def slc_substrate_ontology():

    return _slc_raw(table = 1, sheet = 2)


def slc_chebi_mapping():

    return _slc_raw(table = 2, sheet = 1)
