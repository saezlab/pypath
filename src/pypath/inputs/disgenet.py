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

import collections
import csv

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.utils.mapping as mapping


def disgenet_annotations(dataset = 'curated'):
    """
    Downloads and processes the list of all human disease related proteins
    from DisGeNet.
    Returns dict of dicts.

    @dataset : str
        Name of DisGeNet dataset to be obtained:
        `curated`, `literature`, `befree` or `all`.
    """

    DisGeNetAnnotation = collections.namedtuple(
        'DisGeNetAnnotation',
        [
            'disease',
            'score',
            'dsi',
            'dpi',
            'nof_pmids',
            'nof_snps',
            'source',
        ]
    )

    url = urls.urls['disgenet']['url'] % dataset
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
                DisGeNetAnnotation(
                    disease = rec['diseaseName'],
                    score = float(rec['score']),
                    dsi = float(rec['DSI']) if rec['DSI'] else None,
                    dpi = float(rec['DPI']) if rec['DPI'] else None,
                    nof_pmids = int(rec['NofPmids']),
                    nof_snps = int(rec['NofSnps']),
                    source = tuple(
                        x.strip()
                        for x in rec['source'].split(';')
                    ),
                )
            )

    return data
