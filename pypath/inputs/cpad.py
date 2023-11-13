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

import csv
import collections

import pypath.share.curl as curl
import pypath.share.common as common
import pypath.resources.urls as urls
import pypath.utils.mapping as mapping


def get_cpad():

    url = urls.urls['cpad']['url']
    c = curl.Curl(url, silent = False, large = True, encoding = 'iso-8859-1')
    reader = csv.DictReader(c.result, delimiter = '\t')

    return reader


def cpad_annotations(include_unknown_type = False):

    CpadAnnotation = collections.namedtuple(
        'CpadAnnotation',
        [
            'regulator_type',
            'effect_on_pathway',
            'pathway',
            'effect_on_cancer',
            'effect_on_cancer_outcome',
            'cancer',
            'pathway_category',
        ]
    )

    cpad = get_cpad()

    result = collections.defaultdict(set)

    for rec in cpad:
        if rec['Regulator'] == 'NULL':
            continue

        for regulator in rec['Regulator'].split(' and '):
            uniprot = mapping.map_name0(regulator, 'genesymbol', 'uniprot')

            if uniprot:
                regulator_name = uniprot
                regulator_type = 'protein'

            else:
                mirbase = mapping.map_name(
                    'hsa-%s' % regulator,
                    'mir-mat-name',
                    'mirbase',
                )

                if not mirbase:
                    mirbase = mapping.map_name(
                        'hsa-%s' % regulator,
                        'mir-name',
                        'mirbase',
                    )

                if mirbase:
                    regulator_name = mirbase
                    regulator_type = 'mirna'

                else:
                    if include_unknown_type:
                        regulator_name = regulator
                        regulator_type = 'unknown'

                    else:
                        continue

            if isinstance(regulator_name, str):
                regulator_name = (regulator_name,)

            for regulator_name_0 in regulator_name:
                record = CpadAnnotation(
                    regulator_type = regulator_type,
                    effect_on_pathway = rec['Regulator_Type'],
                    effect_on_cancer = rec['Regulation_Type'],
                    effect_on_cancer_outcome = rec['Outcome_Description'],
                    pathway = rec['Pathway'],
                    pathway_category = rec['Pathway_Category'],
                    cancer = rec['Cancer'],
                )

                result[regulator_name_0].add(record)

    return dict(result)


def cpad_pathway_cancer():
    """
    Collects only the pathway-cancer relationships. Returns sets of records
    grouped in dicts by cancer and by pathway.
    """

    CpadPathwayCancer = collections.namedtuple(
        'CpadPathwayCancer',
        [
            'pathway',
            'cancer',
            'pathway_category',
            'effect_on_cancer',
            'effect_on_cancer_outcome',
        ]
    )

    cpad = get_cpad()

    by_cancer = collections.defaultdict(set)
    by_pathway = collections.defaultdict(set)

    for rec in cpad:

        record = CpadPathwayCancer(
            pathway = rec['Pathway'],
            cancer = rec['Cancer'],
            pathway_category = rec['Pathway_Category'],
            effect_on_cancer = rec['Regulation_Type'],
            effect_on_cancer_outcome = rec['Outcome_Description'],
        )

        by_cancer[rec['Cancer']].add(record)
        by_pathway[rec['Pathway']].add(record)

    return dict(by_cancer), dict(by_pathway)

