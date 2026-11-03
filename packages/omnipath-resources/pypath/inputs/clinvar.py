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

"""
Variant data from the Clinvar database.
"""

from __future__ import annotations

import io
import csv
import ctypes
import collections

import pypath.share.curl as curl
import pypath.resources.urls as urls

csv.field_size_limit(int(ctypes.c_ulong(-1).value // 2))

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
                'review_status',
                'rs',
                'phenotype_ids',
                'phenotypes',
                'otherids',
                'origin',
                'variation_id',
                'assembly',
                'chromosome',
                'chromosome_accession',
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

        phenotype_ids = tuple(row['PhenotypeIDS'].replace('|', ';').replace(',', ';').split(';'))
        phenotypes = tuple(row['PhenotypeList'].replace('|', ';').replace(',', ';').split(';'))
        otherids = tuple(row['OtherIDs'].replace('|', ';').replace(',', ';').split(';'))

        variant = Variant(
            allele = row['AlleleID'],
            type = row['Type'],
            variant = row['Name'],
            entrez = row['GeneID'],
            genesymbol = row['GeneSymbol'],
            clinical_significance = row['ClinicalSignificance'],
            review_status = row['ReviewStatus'],
            rs = row['RS# (dbSNP)'],
            phenotype_ids = phenotype_ids,
            phenotypes = phenotypes,
            otherids = None if otherids[0] == '-' else otherids,
            origin = row['OriginSimple'],
            variation_id = row['VariationID'],
            assembly = row['Assembly'],
            chromosome = row['Chromosome'],
            chromosome_accession = row['ChromosomeAccession'],
        )
        result.add(variant)

    return list(result)


def clinvar_citations() -> list[tuple]:
    """
    Retrieves citation information of variants

    Returns:
        Citations as a list of named tuples.
    """

    Citation = collections.namedtuple(
        'Citation',
        [
            'allele',
            'variation_id',
            'nsv',
            'citation_source',
            'citation_id'
        ],
        defaults=None
    )

    url = urls.urls['clinvar']['url_citations']

    c = curl.Curl(url, large = True, silent = False)

    response = csv.DictReader(
        c.result,
        delimiter = '\t',
    )

    result = set()

    for row in response:
        
        citation = Citation(
            allele = row['#AlleleID'],
            variation_id = row['VariationID'],
            nsv = row['nsv'],
            citation_source = row['citation_source'],
            citation_id = row['citation_id']
        )

        result.add(citation)

    return list(result)
