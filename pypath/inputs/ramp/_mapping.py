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

import pandas as pd

from ._sqlite import ramp_raw

__all__ = ['ramp_mapping']


def ramp_mapping(
        id_type_a: str,
        id_type_b: str,
        return_df: bool = False,
        curies: bool = False,
    ) -> dict[str, set[str]] | pd.DataFrame:
    """
    Retrieve the mapping between two identifiers.

    Args:
        id_type_a:
            The identifier type of the first identifier.
        id_type_b:
            The identifier type of the second identifier.
        return_df:
            Return a pandas DataFrame instead of a dictionary.
        curies:
            Do not remove CURIEs from the identifiers.

    Returns:
        A dictionary with the mapping between the two identifiers.
    """

    query = (
        'SELECT DISTINCT a.sourceId as id_type_a, b.sourceId as id_type_b '
        'FROM '
        '   (SELECT sourceId, rampId '
        '    FROM source '
        f'   WHERE geneOrCompound = "compound" AND IDtype = "{id_type_a}") a '
        'JOIN '
        '   (SELECT sourceId, rampId '
        '    FROM source '
        f'   WHERE geneOrCompound = "compound" AND IDtype = "{id_type_b}") b '
        'ON a.rampId = b.rampId;'
    )

    con = ramp_raw(tables = 'source', sqlite = True)
    df = pd.read_sql_query(query, con)

    if not curies:

        df[df.columns] = df[df.columns].apply(
            lambda y: [x.split(':', maxsplit = 1)[-1] for x in y],
        )

    return (
        df
            if return_df else
        df.groupby('id_type_a')['id_type_b'].apply(set).to_dict()
    )
