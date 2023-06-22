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

from typing import Generator, NamedTuple

import collections

import pandas as pd

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.inputs.common as inputs_common


_notavail = lambda x: None if x == 'Not Available' else x
_synonyms = lambda x: (
    tuple(sorted(y.strip() for y in x.split('|'))) if x else ()
)


class AdrecsAdr(NamedTuple):
    adr_class: str
    badd: str


class AdrecsChildParent(NamedTuple):
    child: AdrecsAdr
    parent: AdrecsAdr


class AdrecsDrugAdr(NamedTuple):
    drug_badd: str
    drug: str
    adr_badd: str
    adr: str


class AdrecsTerm(NamedTuple):
    adrecs_class: str
    badd: str
    name: str
    synonyms: tuple[str]
    meddra: str


class AdrecsDrug(NamedTuple):
    badd: str
    drug: str
    synonyms: str
    drugbank: str
    pubchem_cid: str
    mesh: str
    kegg: str
    tdd: str


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
        record = AdrecsDrug,
        cell_range = 'A1:H2527',
        synonym_idx = [2],
        return_df = return_df,
    )


def adrecs_adr_ontology(return_df: bool = False) -> list[AdrecsTerm] | pd.DataFrame:
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
        record = AdrecsTerm,
        cell_range = 'A1:E13856',
        synonym_idx = [3],
        return_df = return_df,
    )


def _adrecs_base(
        url_key: str,
        record: str | type,
        cell_range: str,
        synonym_idx: list[int],
        fields: tuple[str] | None = None,
        return_df: bool = False,
    ) -> list[tuple] | pd.DataFrame:

    if isinstance(record, str):

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
    ) -> Generator[AdrecsDrugAdr] | pd.DataFrame:
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

    url = urls.urls['adrecs']['adrecs_drugs']
    c = curl.Curl(url, large = True, silent = False)
    _ = next(c.result)

    for line in c.result:

        yield AdrecsDrugAdr(*line.strip().split('\t'))


def adrecs_hierarchy() -> set[AdrecsChildParent]:
    """
    Child-parent relationships between AdReCS ontology terms.

    Return:
        Set of tuples representing child-parent relationship. Both the child
        and parent terms present with their numeric class and BADD identifiers.
    """

    adr_ontology = adrecs_adr_ontology()

    child_adrs = {
        record.adrecs_class: record.badd
        for record in adr_ontology
    }

    result = set()

    for field in adr_ontology:

        if '.' not in field.adrecs_class:

            continue

        parent_adrecs = field.adrecs_class.rsplit('.', 1)[0]

        result.add(
            AdrecsChildParent(
                child = AdrecsAdr(
                    adr_class = field.adrecs_class,
                    badd = field.badd,
                ),
                parent = AdrecsAdr(
                    adr_class = parent_adrecs,
                    badd = child_adrs.get(parent_adrecs),
                ),
            )
        )

    return result
