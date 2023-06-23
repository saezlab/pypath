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

import collections

import pypath.resources.urls as urls
import pypath.share.curl as curl


def compath_mappings(
        source_db: Literal['kegg', 'wikipathways', 'reactome'] | None = None,
        target_db: Literal['kegg', 'wikipathways', 'reactome'] | None = None,
    ) -> Generator[tuple] | pd.DataFrame:
    """
    Cross-database pathway to pathway mappings from Compath.

    Compath contains proposed and accepted mappings by the users/curators
    between pairs of pathways across databases. The source and target
    databases specify the direction of the mapping.

    Args:
        source_db:
            Name of the source database.
        target_db:
            Name of the target database.
        return_df:
            Return a pandas data frame.

    Returns:
        Tuples of pathway-to-pathway mappings.
    """

    result = _compath_mappings(source_db, target_db)

    return pd.DataFrame(result) if return_df else result


def _compath_mappings(
        source_db: Literal['kegg', 'wikipathways', 'reactome'] | None = None,
        target_db: Literal['kegg', 'wikipathways', 'reactome'] | None = None,
    ) -> Generator[tuple]:

    url = urls.urls['compath']['url']
    c = curl.Curl(url, large = True)

    result = set()
    fields = (
        'pathway1',
        'pathway_id_1',
        'source_db',
        'relation',
        'pathway2',
        'pathway_id_2',
        'target_db',
    )
    record = collections.namedtuple('CompathPathwayToPathway', fields)

    for line in c.result:

        line = line.strip().split('\t')

        if (
            source_db is None or l[2] == source_db and
            target_db is None or l[6] == target_db
        ):

            for db_i, pw_i in zip((2, 6), (1, 5)):

                if line[db_i] == 'kegg':

                    line[pw_i] = line[pw_i][5:]

            yield record(*line)

