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

from typing import Generator, Literal

import json
import re
import collections

import pandas as pd

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.share.session as session

_log = session.Logger(name = 'opentargets_input')._log


def opentargets_general(
        dataset: Literal[
            'assoc_direct',
            'assoc_indirect',
            'adr',
            'expr',
            'associationByOverallIndirect',
            'associationByOverallDirect',
            'fda/significantAdverseDrugReactions',
            'baselineExpression',
        ],
        return_df: bool = False,
        by: str | bool = False,
    ) -> Generator[dict] | dict[str, list[dict]] | pd.DataFrame:
    """
    Download data from the Open Targets database.

    Args:
        dataset:
            Name of a dataset, either as a shorthand synonym or as it is shown
            in the URL.
        return_df:
            Return a pandas data frame.
        by:
            Name of the variable to be used as top level key in the returned
            dictionary. If True, the default grouping variable for the given
            dataset will be used. If False, no grouping will be performed.
    """

    by_defaults = {
        'associationByOverallIndirect': 'diseaseId',
        'associationByOverallDirect': 'diseaseId',
        'fda/significantAdverseDrugReactions': 'chembl_id',
        'baselineExpression': 'id',
    }

    by = by_defaults[dataset] if by == True else by

    result = _opentargets_general(dataset)


    if return_df:

        result = pd.DataFrame(result)

    elif by:

        grouped = collections.defaultdict(list)

        for it in result:

            key = it[by]
            del it[by]
            grouped[key].append(it)

        result = grouped

    return result


def _opentargets_general(
        dataset: Literal[
            'assoc_direct',
            'assoc_indirect',
            'adr',
            'expr',
            'associationByOverallIndirect',
            'associationByOverallDirect',
            'fda/significantAdverseDrugReactions',
            'baselineExpression',
        ],
    ) -> Generator[dict]:

    datasets = {
        'assoc_indirect': 'associationByOverallIndirect',
        'assoc_direct': 'associationByOverallDirect',
        'adr': 'fda/significantAdverseDrugReactions',
        'expr': 'baselineExpression',
    }

    dataset = datasets.get(dataset, dataset)


    url = urls.urls['opentargets']['url'] % dataset
    c = curl.Curl(url, silent = False, large = False)

    repart = re.compile(r'"(part.*\.json)"')
    json_files = repart.findall(c.result)
    url += '/%s'

    for json_name in json_files:

        c = curl.Curl(url % json_name, silent = False, large = True)

        for line in c.result:

            if not line:

                continue

            try:

                contents = json.loads(line)

            except json.JSONDecodeError:

                err = f'Failed to parse JSON from Open Targets data:\n{line}'
                _log(err)
                continue

            yield contents


def opentargets_indirect_score(
        return_df: bool = False,
        by: str | bool = False,
    ) -> Generator[dict] | dict | pd.DataFrame:
    """
    Indirect target-disease association scores from Open Targets.

    Args:
        return_df:
            Return a pandas data frame.
        by:
            Name of the variable to be used as top level key in the returned
            dictionary. If True, the default grouping variable for the given
            dataset will be used. If False, no grouping will be performed.

    Returns:
        Target-disease association records as a list of dicts by default; or
        a pandas data frame if `return_df` is True; or a dict of list of dicts
        if by is not False.
    """

    return opentargets_general('assoc_indirect', return_df, by)


def opentargets_direct_score(
        return_df: bool = False,
        by: bool = False,
    ) -> Generator[dict] | dict | pd.DataFrame:
    """
    Direct target-disease association scores from Open Targets.

    Args:
        return_df:
            Return a pandas data frame.
        by:
            Name of the variable to be used as top level key in the returned
            dictionary. If True, the default grouping variable for the given
            dataset will be used. If False, no grouping will be performed.

    Returns:
        Target-disease association records as a list of dicts by default; or
        a pandas data frame if `return_df` is True; or a dict of list of dicts
        if by is not False.
    """

    return opentargets_general('assoc_direct', return_df, by)


def opentargets_adverse_reactions(
        return_df: bool = False,
        by: bool = False,
    ) -> Generator[dict] | dict | pd.DataFrame:
    """
    Drug adverse reactions from Open Targets.

    Args:
        return_df:
            Return a pandas data frame.
        by:
            Name of the variable to be used as top level key in the returned
            dictionary. If True, the default grouping variable for the given
            dataset will be used. If False, no grouping will be performed.

    Returns:
        Drug-adverse reaction records as a list of dicts by default; or
        a pandas data frame if `return_df` is True; or a dict of list of dicts
        if by is not False.
    """

    return opentargets_general('adr', return_df, by)


def opentargets_baseline_expression(
        return_df: bool = False,
        by: bool = False,
    ) -> Generator[dict] | dict | pd.DataFrame:
    """
    Baseline expression from Open Targets.

    Args:
        return_df:
            Return a pandas data frame.
        by:
            Name of the variable to be used as top level key in the returned
            dictionary. If True, the default grouping variable for the given
            dataset will be used. If False, no grouping will be performed.

    Returns:
        Baseline expression records as a list of dicts by default; or
        a pandas data frame if `return_df` is True; or a dict of list of dicts
        if by is not False.
    """

    return opentargets_general('expr', return_df, by)
