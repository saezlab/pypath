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
#           Ömer Kaan Vural
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

from __future__ import annotations

from typing import Literal
import collections

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.common as common
import pypath.inputs.common as inputs_common
import pypath.utils.taxonomy as taxonomy


def oma_orthologs(
        organism_a: str | int,
        organism_b: str | int,
        rel_type: set[Literal['1:1', '1:n', 'm:1', 'm:n']] | None = None,
        score: float = None
    ) -> list[tuple]:
    """
    Retrieves pairwise relations between two genomes from the OMA
    (Orthologous Matrix) database (https://omabrowser.org/oma/home/).

    Args:
        genome_id_a:
            Name or NCBI Taxonomy ID of the first organism.
        genome_id_b:
            Name or NCBI Taxonomy ID of the second organism.
        rel_type:
            Restrict relations to certain types.
        score:
            Resemblance metric.

    Returns:
        A list with tuples of pairwise orthologous relationships.
    """

    OmaOrthologs = collections.namedtuple(
            'OmaOrthologs',
            (
                'id_a',
                'id_b',
                'rel_type',
                'score',
            ),
        )

    organism_a = taxonomy.ensure_ncbi_tax_id(organism_a)
    organism_b = taxonomy.ensure_ncbi_tax_id(organism_b)
    rel_type = common.to_set(rel_type)
    url = urls.urls['oma']['url']
    page = 1
    n_pages = 1e6
    # first decleration is set to prevent recurrency.
    # But at the end it will return as a list
    result = set()

    while True:

        page_url = f'{url}{organism_a}/{organism_b}/?page={page}&per_page=1000'
        c = curl.Curl(page_url, silent = False)

        if not c.result: break

        c.get_headers()
        n_pages = float(c.resp_headers_dict.get('x-total-count', 1e8)) / 100
        page += 1
        data = inputs_common.json_read(c.result)

        result.update({
            OmaOrthologs(
                id_a = i['entry_1']['canonicalid'],
                id_b = i['entry_2']['canonicalid'],
                rel_type = i['rel_type'],
                score = float(i['score']),
            )
            for i in data
            if (
                (not score or i['score'] >= score) and
                (not rel_type or i['rel_type'] in rel_type)
            )
        })

        if page > n_pages: break

    return list(result)
