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
import csv
import collections
import itertools

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.utils.mapping as mapping


def connectomedb_interactions():
    """
    Retrieves ligand-receptor interactions from connectomeDB2020
    https://asrhou.github.io/NATMI/
    """

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
    resemicol = re.compile(r'; ?')

    url = urls.urls['connectomedb2020']['url']
    c = curl.Curl(url, large = True, silent = False)
    tab = list(csv.DictReader(c.result))

    return [
        ConnectomedbInteraction(
            ligand = row['Ligand gene symbol'],
            ligand_location = resemicol.split(row['Ligand location']),
            receptor = row['Receptor gene symbol'],
            references = rea.findall(row['PMID support']),
        )
        for row in tab
    ]


def connectomedb_annotations():
    """
    Retrieves ligand and receptor annotations from connectomeDB2020
    https://asrhou.github.io/NATMI/
    """

    ConnectomedbAnnotation = collections.namedtuple(
        'ConnectomedbAnnotation',
        [
            'role',
            'location',
        ]
    )

    interactions = connectomedb_interactions()

    result = collections.defaultdict(set)

    for ia in interactions:

        for role in ('ligand', 'receptor'):

            loc_attr = '%s_location' % role
            locations = (
                getattr(ia, loc_attr)
                    if hasattr(ia, loc_attr) else
                ('plasma membrane',)
            )

            uniprots = mapping.map_name(
                getattr(ia, role),
                'genesymbol',
                'uniprot',
            )

            for uniprot, location in itertools.product(uniprots, locations):

                result[uniprot].add(
                    ConnectomedbAnnotation(
                        role = role,
                        location = location,
                    )
                )

    return dict(result)
