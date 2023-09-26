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

from typing import Literal
import collections
import itertools

import pandas as pd

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.common as common
import pypath.inputs.common as inputs_common
import pypath.utils.taxonomy as taxonomy
import pypath.utils.mapping as mapping


def oma_orthologs(
        organism_a: str | int = 'human',
        organism_b: str | int = 'mouse',
        id_type: str | None = None,
        rel_type: set[Literal['1:1', '1:n', 'm:1', 'm:n']] | None = None,
        score: float = None,
        return_df: bool = False,
    ) -> list[tuple] | pd.DataFrame:
    """
    Retrieves pairwise relations between two genomes from the OMA
    (Orthologous Matrix) database (https://omabrowser.org/oma/home/).

    Args:
        organism_a:
            Name or NCBI Taxonomy ID of the first organism.
        organism_b:
            Name or NCBI Taxonomy ID of the second organism.
        id_type:
            OMA by default uses UniProt entry IDs and sometimes other
            identifiers for genes. Set this parameter to control which
            ID type all the identifiers are to be mapped to.
        rel_type:
            Restrict relations to certain types.
        score:
            Lower threshold for similarity metric.
        return_df:
            If True, returns a data frame instead of a list of tuples.

    Returns:
        A list with tuples of pairwise orthologous relationships or a data
        frame with the same records.
    """

    OmaGene = collections.namedtuple(
        'OmaGene',
        (
            'id',
            'oma_group',
            'hog',
            'taxon',
            'chr',
            'start',
            'end',
            'strand',
            'main_isoform',
        )
    )

    OmaOrthology = collections.namedtuple(
            'OmaOrthology',
            (
                'a',
                'b',
                'rel_type',
                'dist',
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

        for rec in data:

            if (
                (score and rec['score'] < score) or
                (rel_type and rec['rel_type'] not in rel_type)
            ):
                continue

            a, b = (
                [
                    OmaGene(
                        id = id_,
                        oma_group = e['oma_group'],
                        hog = e['oma_hog_id'],
                        taxon = e['species']['taxon_id'],
                        chr = e['chromosome'],
                        start = int(e['locus']['start']),
                        end = int(e['locus']['end']),
                        strand = int(e['locus']['strand']),
                        main_isoform = e['is_main_isoform'],

                    )
                    for id_ in _id_translate(
                        id_ = e['canonicalid'],
                        taxon = e['species']['taxon_id'],
                        id_type = id_type,
                    )
                ]
                for e in (rec[f'entry_{ei}'] for ei in (1, 2))
            )


            result.update(
                {
                    OmaOrthology(
                        a = _a,
                        b = _b,
                        rel_type = rec['rel_type'],
                        dist = float(rec['distance']),
                        score = float(rec['score']),
                    )
                    for _a, _b in itertools.product(a, b)
                }
            )

        if page > n_pages: break

    result = list(result)

    if return_df:

        result = pd.DataFrame(result)
        result = pd.concat(
            [
                pd.DataFrame(result.a.tolist()).add_suffix('_a', axis = 1),
                pd.DataFrame(result.b.tolist()).add_suffix('_b', axis = 1),
                result.drop(['a', 'b'], axis = 1),
            ],
            axis = 1,
        )

    return result


def oma_table(
        organism_a: str | int = 'human',
        organism_b: str | int = 'mouse',
        id_type: str | None = None,
        rel_type: set[Literal['1:1', '1:n', 'm:1', 'm:n']] | None = None,
        score: float = None,
        return_df: bool = False,
    ) -> dict[str[set[str]]] | pd.DataFrame:
    """
    Translation table of orthologous gene pairs between two organisms from the
    OMA (Orthologous Matrix) database (https://omabrowser.org/oma/home/).

    Args:
        organism_a:
            Name or NCBI Taxonomy ID of the first organism.
        organism_b:
            Name or NCBI Taxonomy ID of the second organism.
        id_type:
            OMA by default uses UniProt entry IDs and sometimes other
            identifiers for genes. Set this parameter to control which
            ID type all the identifiers are to be mapped to.
        rel_type:
            Restrict relations to certain types.
        score:
            Lower threshold for similarity metric.
        return_df:
            If True, returns a data frame instead of a list of tuples.

    Returns:
        A dict with source organism identifiers as keys and sets of target
        organism identifiers as values; or a two column data frame with
        source organism-target organism identifier pairs.
    """

    full = oma_orthologs(**locals())

    if return_df:

        result = full[['id_a', 'id_b']]

    else:

        result = collections.defaultdict(set)

        for o in full:

            result[o.a.id].add(o.b.id)

    return result


def _id_translate(id_: str, taxon: int, id_type: str | None) -> set[str]:

    if not id_type: return {id_}

    s_id_type = (
        'ensg'
            if id_.startswith('ENS') else
        'uniprot-entry'
            if '_' in id_ else
        'uniprot'
    )

    uniprots = mapping.map_name(
        id_,
        s_id_type,
        'uniprot',
        ncbi_tax_id = taxon,
    )

    return mapping.map_names(
        uniprots,
        'uniprot',
        id_type,
        ncbi_tax_id = taxon,
    ) if uniprots else set()
