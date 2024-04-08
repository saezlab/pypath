#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright 2014-2024
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

import json

import pandas as pd

import pypath.resources.urls as urls
import pypath.share.curl as curl


def idtypes(
        pairs: bool = True,
        raw: bool = False,
    ) -> dict[str, pd.DataFrame] | set[tuple[str, str]] | dict:
    """
    Identifier types in the UniProt ID mapping service.

    Args:
        pairs:
            Process the data into pairs of identifiers.
        raw:
            Return the raw data as extracted from JSON.

    Returns:
        The JSON contents as a dict if `raw` is `True`,
        a list of tuples if `pairs` is `True`,
        otherwise a set of tuples of ID types.
    """

    url = urls.urls['uniprot_idmapping']['fields']
    c = curl.Curl(url, large = False, silent = False)
    data = json.loads(c.result)

    if raw:

        return data

    groups = (
        pd.DataFrame(data['groups']).
        explode('items').
        reset_index(drop = True)
    )
    groups = (
        pd.concat(
            [
                groups['groupName'],
                pd.DataFrame(groups['items'].tolist())
            ],
            axis = 1,
        ).
        rename(columns = {'from': 'from_'})
    )

    rules = pd.DataFrame(data['rules'])

    if not pairs:

        return {'groups': groups, 'rules': rules}

    rules = {int(r.ruleId): r.tos for r in rules.itertuples()}
    groups.fillna(-1., inplace = True)

    result = set()

    for idtype in groups.itertuples():

        tos = rules.get(int(idtype.ruleId), [])

        from_to = {(idtype.name, t) for t in tos}

        if idtype.from_:

            result.update(from_to)

        if idtype.to:

            result.update({t[::-1] for t in from_to})

    return result
