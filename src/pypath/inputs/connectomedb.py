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

import re
import csv
import collections

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.utils.mapping as mapping


def connectomedb_interactions():

    ConnectomedbInteraction = collections.namedtuple(
        'ConnectomedbInteraction',
        [
            'ligand',
            'ligand_location',
            'receptor',
            'references',
        ]
    )

    rea = re.compile(r'<a[^>]+>([^<]*)</a>')

    url = urls.urls['connectomedb2020']['url']
    c = curl.Curl(url, large = True, silent = False)
    tab = list(csv.DictReader(c.result))

    return [
        ConnectomedbInteraction(
            ligand = row['Ligand gene symbol'],
            ligand_location = row['Ligand location'],
            receptor = row['Receptor gene symbol'],
            references = rea.findall(row['PMID support']),
        )
        for row in tab
    ]



