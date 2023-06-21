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

import bs4

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.utils.mapping as mapping


def dgidb_interactions() -> list[tuple]:
    """
    Retrieves drug-gene interactions from DGIdb.

    Returns:
        A list with tuples. Tuples are dgidb interactons
    """

    result = set()

    DgidbInteraction = collections.namedtuple(
        'DgidbInteraction',
        [
            'genesymbol',
            'entrez',
            'resource',
            'type',
            'drug_name',
            'drug_chembl',
            'score',
            'pmid'
        ],
    )

    url = urls.urls['dgidb']['interactions']
    c = curl.Curl(url = url, silent = False, large = True)
    interactions = csv.DictReader(c.result, delimiter = '\t')

    for interaction in interactions:

        interaction = {k: v or None for k, v in interaction.items()}

        dgidb_interaction = DgidbInteraction(
            genesymbol = interaction['gene_name'],
            entrez = interaction['entrez_id'],
            resource = interaction['interaction_claim_source'],
            type = interaction['interaction_types'],
            drug_name = interaction['drug_claim_primary_name'],
            drug_chembl = interaction['drug_concept_id'],
            score = interaction['interaction_group_score'],
            pmid = interaction['PMIDs'],
        )

        result.add(dgidb_interaction)

    return list(result)


def dgidb_annotations():
    """
    Downloads druggable protein annotations from DGIdb.
    """

    DgidbAnnotation = collections.namedtuple(
        'DgidbAnnotation',
        ['category'],
    )


    url = urls.urls['dgidb']['categories']
    c = curl.Curl(url = url, silent = False, large = True)
    data = csv.DictReader(c.result, delimiter = '\t')

    result = collections.defaultdict(set)

    for rec in data:

        uniprots = mapping.map_name(
            rec['entrez_gene_symbol'],
            'genesymbol',
            'uniprot',
        )

        for uniprot in uniprots:
            result[uniprot].add(
                DgidbAnnotation(
                    category = rec['category']
                )
            )

    return dict(result)


def get_dgidb_old():
    """
    Deprecated. Will be removed soon.

    Downloads and processes the list of all human druggable proteins.
    Returns a list of GeneSymbols.
    """

    genesymbols = []
    url = urls.urls['dgidb']['main_url']
    c = curl.Curl(url, silent = False)
    html = c.result
    soup = bs4.BeautifulSoup(html, 'html.parser')
    cats = [
        o.attrs['value']
        for o in soup.find('select', {'id': 'gene_categories'})
        .find_all('option')
    ]

    for cat in cats:
        url = urls.urls['dgidb']['url'] % cat
        c = curl.Curl(url)
        html = c.result
        soup = bs4.BeautifulSoup(html, 'html.parser')
        trs = soup.find('tbody').find_all('tr')
        genesymbols.extend([tr.find('td').text.strip() for tr in trs])

    return mapping.map_names(genesymbols, 'genesymbol', 'uniprot')
