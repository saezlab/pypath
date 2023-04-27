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

import sys
sys.path.append('/home/exdeval/.stuff/biolab/star/pypath/pypath')

from pypath.share import curl
from collections import namedtuple
from resources.urls import urls

to_numeric = {
    'z_score',
    'confidence_score',
    'confidence_score'
}

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
        
        for key, value in data.items():
            if key in to_numeric:
                data[key] = str_to_num(value)
            elif key == 'source_score':
                new_value = value.split('=')[1]
                new_value = str_to_num(new_value)
                data[key] = new_value

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


def str_to_num(string):
    try:
        return int(string)
    except ValueError:
        return float(string)