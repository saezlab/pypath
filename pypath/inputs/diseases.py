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
#           Sebastian Lobentanzer
#           Erva Ulusoy
#           Olga Ivanova
#           Ahmet Rifaioglu
#           Melih Darcan
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

from pypath.share import curl
from collections import namedtuple
from pypath.resources.urls import urls


def diseases_general(query, filtered=False):

    query_type = 'filtered' if filtered else 'full'

    url = urls['diseases']['url'] % (query, query_type)

    fieldnames = [
        'gene_identifier',
        'gene_name',
        'disease_identifier',
        'disease_name',
    ]

    if query == 'textmining':
        fieldnames.extend(
            [
                'z_score',
                'confidence_score',
                'url'
            ]
        )

    elif query == 'knowledge':
        fieldnames.extend(
            [
                'source_database',
                'evidence_type',
                'confidence_score'
            ]
        )

    elif query == 'experiments':
        fieldnames.extend(
            [
                'source_database',
                'source_score',
                'confidence_score'
            ]
        )

    elif query == 'integrated':
        print('Not supported yet.')
        exit(0)

    else:
        print('Problem in function call. Check arguments.')
        exit(1)

    c = curl.Curl(
        url,
        silent=False,
        large=True,
        encoding="utf-8",
        default_mode="r",
    )

    Interaction = namedtuple('DISEASESInteraction', fieldnames)
    interactions = list()

    for line in c.result:

        line = line.strip('\n ')
        data = line.split("\t")

        if data[-1] == "\n":
            del data[-1]

        data = {
            fieldname: element if element != "" else None
            for (fieldname, element) in zip(fieldnames, data)
        }

        interactions.append(Interaction(**data))

    return interactions


def textmining_full():
    return diseases_general('textmining', filtered=False)


def textmining_filtered():
    return diseases_general('textmining', filtered=True)


def knowledge_full():
    return diseases_general('knowledge', filtered=False)


def knowledge_filtered():
    return diseases_general('knowledge', filtered=True)


def experiments_full():
    return diseases_general('experiments', filtered=False)


def experiments_filtered():
    return diseases_general('experiments', filtered=True)
