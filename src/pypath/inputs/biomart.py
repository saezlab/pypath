#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#  Helps to translate from the mouse data to human data
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

import re
import os

import pypath.share.common as common
import pypath.share.curl as curl
import pypath.resources.urls as urls


def biomart_query(attr, transcript = False):

    # attr = 'hgnc_symbol'

    rewsp = re.compile(r'\n\s+')

    xml_template_path = os.path.join(common.DATA, 'ensembl_biomart_query.xml')

    with open(xml_template_path, 'r') as fp:

        xml_template = fp.read()

    ens_id_type = 'transcript' if transcript else 'gene'

    xml_query = xml_template % (
        ens_id_type,
        attr,
    )
    xml_query = rewsp.sub('', xml_query)

    biomart_url = urls.urls['ensembl']['biomart_url'] % xml_query

    c = curl.Curl(biomart_url, large = True, silent = False)

    for line in c.result:

        line = line.strip('\n\r').split('\t')

        if len(line) >= 2:

            yield line