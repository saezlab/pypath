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

import re
import collections

import bs4

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.utils.taxonomy as taxonomy


def lncrnadb_interactions():

    LncrnadbInteraction = collections.namedtuple(
        'LncrnadbInteraction',
        (
            'lncrna',
            'partner',
            'type',
            'organism',
            'pmid',
        ),
    )


    renondigit = re.compile(r'[^\d]+')

    url = urls.urls['lncrnadb']['url_rescued']
    c = curl.Curl(
        url,
        silent = False,
        large = True,
        encoding = 'utf-8',
    )

    b = bs4.BeautifulSoup(c.fileobj, features = 'xml')

    result = []

    for res in b.findAll('results'):

        lncrna = res.find('nomenclature').find('name').text

        for sp in res.find('species').findAll('entry'):

            organism = sp.attrs['species'].split('(')[0].strip()

            for assoc in res.find('association').findAll('association'):

                partner  = assoc.find('componentid').text
                typ      = assoc.find('componenttype').text.lower()
                pmid     = renondigit.sub('', assoc.find('pubmedid').text)

                result.append(
                    LncrnadbInteraction(
                        lncrna = lncrna,
                        partner = partner,
                        type = typ,
                        organism = taxonomy.ensure_ncbi_tax_id(organism),
                        pmid = pmid,
                    )
                )

    return result
