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
#           Tennur Kılıç
#           Ömer Kaan Vural
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

"""
Variant data from the Clinvar database.
"""

from __future__ import annotations

import io
import csv
import sys
import collections

import pypath.share.curl as curl
import pypath.resources.urls as urls

csv.field_size_limit(sys.maxsize)


def clinvar_raw() -> list[tuple]:
    """
    Retrieves variant data from the Clinvar database.

    Returns:
        Variants as a list of named tuples.
    """

    Variant = collections.namedtuple(
            'Variant',
            [
                'allele',
                'type',
                'variant',
                'entrez',
                'genesymbol',
                'clinical_significance',
                'rs',
                'phenotype_ids',
                'phenotypes',
                'origin',
                'variation_id',
            ],
            defaults = None
        )

    url = urls.urls['clinvar']['url']
    c = curl.Curl(url, large = True, silent = False)
    c.gzfile.seek(1) # get rid of a stray `#` character

    response = csv.DictReader(
        io.TextIOWrapper(c.gzfile),
        dialect = 'excel-tab',
    )

    result = set()

    for row in response:

        phenotype_ids = tuple(row['PhenotypeIDS'].replace('|', ';').split(';'))
        phenotypes = tuple(row['PhenotypeList'].replace('|', ';').split(';'))

        variant = Variant(
            allele = row['AlleleID'],
            type = row['Type'],
            variant = row['Name'],
            entrez = row['GeneID'],
            genesymbol = row['GeneSymbol'],
            clinical_significance = row['ClinicalSignificance'],
            rs = row['RS# (dbSNP)'],
            phenotype_ids = phenotype_ids,
            phenotypes = phenotypes,
            origin = row['OriginSimple'],
            variation_id = row['VariationID']
        )
        result.add(variant)

    return list(result)
