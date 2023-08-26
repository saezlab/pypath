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

from collections import namedtuple

import json
import bs4
import requests
import json

import pandas as pd

import pypath.share.curl as curl
import pypath.resources.urls as urls


def ddinter_n_drugs() -> int:
    """
    Retrieves number of drug records in the DDInter database.
    """

    cookies = {
        'csrftoken': '975RT1K1FXpLX6wZdcyTMIRdOAMwNNkx1MNFAbeJrlrRBFI5DbNmaQPz7tj2UNFZ',
    }

    data = {
        'csrfmiddlewaretoken': 'ToKztkgeCHgBfAbezQL7GULrQABtRmHgL3snauKWo5iHT9nkZP0A42JN9t8ZYm2I',
        'draw': '1',
        'columns[0][data]': '',
        'columns[0][name]': '',
        'columns[0][searchable]': 'true',
        'columns[0][orderable]': 'false',
        'columns[0][search][value]': '',
        'columns[0][search][regex]': 'false',
        'start': '0',
        'length': '24',
        'search[value]': '',
        'search[regex]': 'false',
        'custom-length': '25',
    }

    url = urls.urls['ddinter']['source']
    response = requests.post(
        url,
        cookies = cookies,
        data = data,
        verify = False,
    )

    return int(response.json()['recordsTotal'])


def ddinter_identifiers(drug: str) -> list[tuple]:
    """
    DrugBank, ChEMBL and PubChem identifiers of a drug in DDInter.

    Args:
        drug:
            DDInter drug identifier.

    Returns:
        list of mapping namedtuple (single)
    """
    url = urls.urls['ddinter']['mapping'] % drug

    c = curl.Curl(url, silent = False, large = True)

    soup = bs4.BeautifulSoup(c.fileobj, 'html.parser')
    refs = soup.find_all('a')

    links = [link.get('href', '') for link in refs]
    result = set()

    fields = ('ddinter', 'drugbank', 'chembl', 'pubchem')
    record = namedtuple('DdinterIdentifiers', fields, defaults = (None, ) * len(fields))
    mapping_targets = ['drugbank', 'chembl', 'pubchem']
    mapping_dict = {}

    for link in links:

        for target in mapping_targets:

            if target in link:

                mapping_dict[target] = link.split('/')[-1]
                break

    result.add(record(ddinter = drug, **mapping_dict))

    return list(result)


def ddinter_mappings(return_df: bool = False) -> list[tuple] | pd.DataFrame:
    """
    DrugBank, ChEMBL and PubChem identifiers of all drugs in DDInter.

    Args:
        return_df:
            Return a pandas data frame.

    Returns:
        list of mapping namedtuples
    """

    result = []

    for idx in range(1, ddinter_n_drugs() + 1):

        ddinter_drug = f'DDInter{idx}'
        
        result.extend(ddinter_identifiers(ddinter_drug))

    return pd.DataFrame(result) if return_df else result


def _ensure_hashable(data):

    if isinstance(data, (dict, list, set)):

        return tuple(data)

    return data


def ddinter_drug_interactions(
        drug: str,
        return_df: bool = False,
    ) -> list[tuple] | pd.DataFrame:
    """
    Interactions of one single drug from the DDInter database.

    Args:
        drug:
            A DDInter drug identifier.
        return_df:
            Return a pandas data frame.

    Returns:
        Drug-drug interaction tuples with drug ids, drug names, interaction
        level and actions.
    """

    url = urls.urls['ddinter']['interaction'] % drug
    c = curl.Curl(url)
    data = json.loads(c.result)

    result = set()
    record = namedtuple(
        'DdinterInteraction',
        (
            'drug1_id',
            'drug1_name',
            'drug2_id',
            'drug2_name',
            'level',
            'actions',
        ),
        defaults = None,
    )

    drug1_id = data['info']['id']
    drug1_name = data['info']['Name']

    for interaction_fe in data['interactions']:

        interaction = record(
            drug1_id = drug1_id,
            drug1_name = drug1_name,
            drug2_id = _ensure_hashable(interaction_fe['id']),
            drug2_name =  _ensure_hashable(interaction_fe['name']),
            level = _ensure_hashable(interaction_fe['level']),
            actions = _ensure_hashable(interaction_fe['actions']),
        )

        result.add(interaction)

    return pd.DataFrame(result) if return_df else result


def ddinter_interactions(return_df: bool = False) -> list[tuple] | pd.DataFrame:
    """
    Drug-drug interactions from the DDInter database.

    Args:
        return_df:
            Return a pandas data frame.

    Returns:
        list of interaction namedtuple with drug ids, drug names, interaction level and actions
    """
    result = [
        ddinter_drug_interactions(drug = f'DDInter{i}')
        for i in range(1, ddinter_n_drugs()+1)
    ]

    return pd.DataFrame(result) if return_df else result
