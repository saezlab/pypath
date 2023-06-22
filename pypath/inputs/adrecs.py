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

from typing import Generator

import collections

import pandas as pd

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.inputs.common as inputs_common


_notavail = lambda x: None if x == 'Not Available' else x
_synonyms = lambda x: (
    tuple(sorted(y.strip() for y in x.split('|'))) if x else ()
)


def adrecs_drug_identifiers(
        return_df: bool = False,
    ) -> list[tuple] | pd.DataFrame:
    """
    Drug identifiers from the AdReCS database.

    IUPAC name, synonyms, DrugBank, MeSH, KEGG and TDD IDs of drugs.
    http://www.bio-add.org/ADReCS/index.jsp

    Args:
        return_df:
            Return a pandas data frame.

    Returns:
        List of tuples or data frame of drug identifiers.
    """

    return _adrecs_base(
        url_key = 'drug_information',
        record_name = 'Drug',
        cell_range = 'A1:H2527',
        fields = (
            'badd',
            'drug',
            'synonyms',
            'drugbank',
            'pubchem_cid',
            'mesh',
            'kegg',
            'tdd',
        ),
        synonym_idx = [2],
        return_df = return_df,

    )


def adrecs_adr_ontology(return_df: bool = False) -> list[tuple] | pd.DataFrame:
    """
    Adverse drug reaction (ADR) ontology from the AdReCS database.

    Args:
        return_df:
            Return a pandas data frame.

    Returns:
        List of tuples or data frame of adverse drug reaction terms.
    """

    return _adrecs_base(
        url_key = 'terminology',
        record_name = 'Term',
        cell_range = 'A1:E13856',
        fields = ('adrecs_class', 'badd', 'name', 'synonyms', 'meddra'),
        synonym_idx = [3],
        return_df = return_df,
    )


def _adrecs_base(
        url_key: str,
        record_name: str,
        cell_range: str,
        fields: tuple[str],
        synonym_idx: list[int],
        return_df: bool = False,
    ) -> list[tuple] | pd.DataFrame:

    record = collections.namedtuple(f'Adrecs{record_name}', fields)

    url = urls.urls['adrecs'][url_key]
    path = curl.Curl(url, silent = False, large = True)
    contents = inputs_common.read_xls(path.outfile, cell_range = cell_range)
    result = []

    for line in contents[1:]:

        line = [_notavail(x) for x in line]

        for isyn in synonym_idx:

            line[isyn] = _synonyms(line[isyn])

        result.append(record(*line))

    return pd.DataFrame(result) if return_df else result


def adrecs_drug_adr(
        return_df: bool = False,
    ) -> Generator[tuple] | pd.DataFrame:
    """
    Drug-ADR pairs from the AdReCS database.

    Args:
        return_df:
            Return a pandas data frame.

    Returns:
        List of tuples or data frame of drug-ADR pairs.
    """

    result = _adrecs_drug_adr()

    return pd.DataFrame(result) if return_df else result


def _adrecs_drug_adr():

    record = collections.namedtuple(
        'AdrecsDrugAdr',
        ('drug_badd',  'drug', 'adr_badd', 'adr'),
    )

    url = urls.urls['adrecs']['adrecs_drugs']
    c = curl.Curl(url, large = True, silent = False)
    _ = next(c.result)

    for line in c.result:

        yield record(*line.strip().split('\t'))
        
def adrecs_extract_child_parent_relationship() -> list[tuple]:
    data = adrecs_adr_ontology()
    
    fields = [
        'child_adrecs_id',
        'child_adr_id',
        'parent_adrecs_id',
        'parent_adr_id'
    ]  

    child_adr_ids = {record.adrecs_class: record.badd for record in data}
    
    result = set()
    record = collections.namedtuple('AdrecsChildParentRelationship', fields)
    
    for field in data:
        if '.' not in field.adrecs_class:
            continue
        
        child_adrecs_id = field.adrecs_class
        child_adr_id = field.badd
        parent_adrecs_id = child_adrecs_id.rsplit('.', 1)[0]
        parent_adr_id = child_adr_ids.get(parent_adrecs_id)

        relationship = (
            child_adrecs_id,
            child_adr_id,
            parent_adrecs_id,
            parent_adr_id
        )
        
        result.add(
            record(**dict(zip(fields, relationship)))
        )
        
    return list(result)


