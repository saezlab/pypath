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

import csv
import collections

import pypath.share.curl as curl
import pypath.share.common as common
import pypath.resources.urls as urls
import pypath.utils.mapping as mapping


def adhesome_interactions():

    AdhesomeInteraction = collections.namedtuple(
        'AdhesomeInteraction',
        ['source', 'target', 'effect', 'type', 'pmid'],
    )

    url = urls.urls['adhesome']['interactions']

    c = curl.Curl(url, large = True, silent = False)

    data = csv.DictReader(c.result, delimiter = ',')

    result = []

    for rec in data:

        result.append(
            AdhesomeInteraction(
                source = rec['Source'],
                target = rec['Target'],
                effect = rec['Effect'],
                type   = common.upper0(rec['Type']),
                pmid   = rec['PMID'],
            )
        )

    return result


def adhesome_annotations():

    AdhesomeAnnotation = collections.namedtuple(
        'AdhesomeAnnotation',
        ['mainclass', 'intrinsic'],
    )

    result = collections.defaultdict(set)

    url = urls.urls['adhesome']['components']
    c = curl.Curl(url, large = True, silent = False)

    data = csv.DictReader(c.result, delimiter = ',')

    for rec in data:
        uniprots = rec['Swiss-Prot ID']

        for uniprot in uniprots.split(','):
            uniprot = uniprot.strip()

            if uniprot == 'null':
                continue

            for _uniprot in mapping.map_name(uniprot, 'uniprot', 'uniprot'):
                result[uniprot].add(AdhesomeAnnotation(
                    mainclass = (
                        common.upper0(rec['Functional Category'].strip())
                    ),
                    intrinsic = rec['FA'].strip() == 'Intrinsic Proteins',
                ))

    return result
