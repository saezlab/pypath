#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2023
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: Dénes Türei (turei.denes@gmail.com)
#           Nicolàs Palacio
#           Sebastian Lobentanzer
#           Erva Ulusoy
#           Olga Ivanova
#           Ahmet Rifaioglu
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

from __future__ import annotations

from typing import Generator, Literal

import collections

import pandas as pd

import pypath.share.curl as curl
import pypath.resources.urls as urls


_molecule_fields =(
    'pathway_id',
    'url',
    'pathway_name',
    'evidence_code',
    'organism'
)

_FIELDS = {
    'uniprot2all': ('uniprot_id',) + _molecule_fields,
    'chebi2all': ('chebi_id',) + _molecule_fields,
    'pathways': (
        'pathway_id',
        'pathway_name',
        'organism',
    ),
    'pathway_relations': (
        'parent',
        'child',
    ),
}

_RECORD_NAMES = {
    'uniprot2all': 'Uniprot',
    'chebi2all': 'Chebi',
    'pathways': 'Pathway',
    'pathway_relations': 'PathwayRelation',
}


def _reactome_data(
        dataset: Literal[
            'uniprot2all',
            'chebi2all',
            'pathways',
            'pathway_relations',
        ],
        return_df: bool = False,
    ) -> Generator[tuple] | pd.DataFrame:
    """
    Retrieves pathway information by given input
    Args:
        source:
            - uniprot2all (fields: uniprot_id, pathway_id, url,
              pathway_name, evidence_code, organism)
            - chebi2all (fields: chebi_id, pathway_id, url,
              pathway_name, evidence_code, organism)
            - pathways (fields: pathway_id, pathway_name, organism)
            - pathway_relations (fields: parent, child)
        return_df:
            Return a pandas data frame.

    Returns:
        Records of the requested dataset.
    """

    result = _reactome_data_gen(dataset)

    return pd.DataFrame(result) if return_df else result


def _reactome_data_gen(
        dataset: Literal[
            'uniprot2all',
            'chebi2all',
            'pathways',
            'pathway_relations',
        ]
    ) -> Generator[tuple]:

    url = urls.urls['reactome'][dataset]
    c = curl.Curl(url, large = True)


    result = set()
    record = collections.namedtuple('Reactome', fields)

    for line in c.result:

        line = line.strip('\n').split('\t')
        yield record(*line)
