#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2022
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: Dénes Türei (turei.denes@gmail.com)
#           Nicolàs Palacio
#           Olga Ivanova
#           Sebastian Lobentanzer
#           Ahmet Rifaioglu
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

import re
import csv
import collections


import pypath.inputs.common as inputs_common
import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.settings as settings
import pypath.utils.mapping as mapping


def panglaodb_raw():
    """
    Cell type marker genes from PanglaoDB (https://panglaodb.se/index.html).
    """

    url = urls.urls['panglaodb']['url']
    headers = [settings.get('user_agent')]
    c = curl.Curl(url, silent = False, large = True, req_headers = headers)

    result = list(csv.DictReader(c.result, delimiter = '\t'))

    return result
