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

from __future__ import annotations

import csv
import collections

import pypath.share.curl as curl
import pypath.inputs.disgenet._auth as _auth
import pypath.inputs.disgenet._records as _records


def disgenet_variants() -> dict[str, _records.VariantGeneMapping]:
    """
    Downloads and processes variant-gene mappings.
    Returns a dict where the \'snpId\' is the key.
    """

    url = urls.urls['disgenet']['variant_gene_mappings']
    c = curl.Curl(
        url,
        silent = False,
        large = True,
        encoding = 'utf-8',
        default_mode = 'r',
    )
    reader = csv.DictReader(c.result, delimiter = '\t')
    mapping = dict()

    for rec in reader:

        snpId = rec.pop('snpId')

        try:
            match = False

            for index, entry in enumerate(mapping[snpId]):

                if (
                    rec['geneId'] == entry['geneId']
                    and rec['geneSymbol'] == entry['geneSymbol']
                ):

                    match = True

                    if isinstance(mapping[snpId][index]['sourceId'], list):

                        mapping[snpId][index]['sourceId'].append(rec['sourceId'])

                    else:

                        mapping[snpId][index]['sourceId'] = [
                            mapping[snpId][index]['sourceId'],
                            rec['sourceId'],
                        ]

                    break

            if not match:

                mapping[snpId].append(rec)

        except KeyError:

            mapping[snpId] = [rec]

    for key, values in mapping.items():

        for index, value in enumerate(values):

            mapping[key][index] = _records.VariantGeneMapping(
                value['geneId'],
                value['geneSymbol'],
                tuple(value['sourceId']),
            )

    return mapping


def disgenet_diseases() -> dict[str, _records.DiseaseIdMapping]:
    """
    Downloads and processes disease-id mappings.
    Returns a dict where the \'diseaseId\' is the key.
    """

    url = urls.urls['disgenet']['disease_id_mappings']
    c = curl.Curl(
        url,
        silent = False,
        large = True,
        encoding = 'utf-8',
        default_mode = 'r',
    )
    reader = csv.DictReader(c.result, delimiter = '\t')
    mapping = dict()

    for rec in reader:

        diseaseId = rec.pop('diseaseId')
        name = rec.pop('name')
        rec = _records.IdType(
            rec['vocabulary'],
            rec['code'],
            rec['vocabularyName'],
        )

        try:

            mapping[diseaseId]['vocabularies'].append(rec)

        except KeyError:

            mapping[diseaseId] = dict()
            mapping[diseaseId]['name'] = name
            mapping[diseaseId]['vocabularies'] = [rec]

    for key, value in mapping.items():

        mapping[key] = _records.DiseaseIdMapping(
            value['name'],
            tuple(value['vocabularies']),
        )

    return mapping


def disgenet_annotations(
        dataset = 'curated',
    ) -> dict[str, _records.DisgenetAnnotation]:
    """
    Downloads and processes the list of all human disease related proteins
    from DisGeNet.
    Returns dict of dicts.

    Args:
        dataset:
            Name of DisGeNet dataset to be obtained:
            `curated`, `literature`, `befree` or `all`.
    """

    url = urls.urls['disgenet']['annotations'] % dataset
    c = curl.Curl(
        url,
        silent = False,
        large = True,
        encoding = 'utf-8',
        default_mode = 'r',
    )
    reader = csv.DictReader(c.result, delimiter = '\t')
    data = collections.defaultdict(set)

    for rec in reader:

        uniprots = mapping.map_name(
            rec['geneSymbol'],
            'genesymbol',
            'uniprot',
        )

        if not uniprots:
            continue

        for uniprot in uniprots:

            data[uniprot].add(
                _records.DisgenetAnnotation(
                    disease = rec['diseaseName'],
                    type = rec['diseaseType'],
                    score = float(rec['score']),
                    dsi = float(rec['DSI']) if rec['DSI'] else None,
                    dpi = float(rec['DPI']) if rec['DPI'] else None,
                    nof_pmids = int(rec['NofPmids']),
                    nof_snps = int(rec['NofSnps']),
                    source = tuple(x.strip() for x in rec['source'].split(';')),
                )
            )

    return dict(data)
