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

from collections import namedtuple

import pandas as pd

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.share.common as common

_NUMERIC_FIELDS = {
    'z_score',
    'confidence',
    'confidence',
    'source_score',
}


def diseases_general(
        data_origin: Literal['textmining',  'knowledge', 'experiments'],
        filtered: bool = False,
        return_df: bool = False,
    ) -> Generator[tuple] | pd.DataFrame:
    """
    Retrieve a dataset from the DISEASES database from Jensen Lab.

    Warning: The "textmining" datasets are enormous!

    Args:
        data_origin:
            The data collection method.
        filtered:
            Download the filtered dataset instead of the full.
        return_df:
            Return a pandas data frame.
    """

    result = _diseases_general(data_origin, filtered)

    return pd.DataFrame(result) if return_df else result


def _diseases_general(
        data_origin: Literal['textmining',  'knowledge', 'experiments'],
        filtered: bool = False,
    ) -> Generator[tuple]:
    """
    Args:
        data_origin:
            The data collection method.
        filtered:
            Download the filtered dataset instead of the full.
    """

    query_type = 'filtered' if filtered else 'full'

    url = urls.urls['diseases']['url'] % (data_origin, query_type)

    query_fields = {
        'textmining':
            [
                'z_score',
                'confidence',
                'url'
            ],
        'knowledge':
            [
                'resource',
                'evidence_type',
                'confidence'
            ],
        'experiments':
            [
                'resource',
                'source_score',
                'confidence'
            ],
    }

    fields = [
        'gene_id',
        'genesymbol',
        'disease_id',
        'disease',
    ] + query_fields[data_origin]


    record = namedtuple('DiseasesInteraction', fields)

    c = curl.Curl(url, silent = False, large = True)
    interactions = list()

    def proc_field(value, key):

        if key == 'source_score':

            value = value.split('=')[1]

        if key in _NUMERIC_FIELDS:
            
            if '.' in value:
                num_type = float
            else:
                num_type = int if common.is_int(value) else float

            value = num_type(value)

        return (key, value)


    for line in c.result:

        line = line.strip('\n ').split('\t')

        line = dict(
            proc_field(value, key)
            for key, value in zip(fields, line)
        )

        yield record(**line)


def textmining_full(
        return_df: bool = False,
    ) -> Generator[tuple] | pd.DataFrame:
    """
    Full textmining dataset of the DISEASES database from Jensen Lab.

    Warning: The "textmining" datasets are enormous!

    Args:
        return_df:
            Return a pandas data frame.
    """

    return diseases_general(
        data_origin = 'textmining',
        filtered = False,
        return_df = return_df,
    )


def textmining_filtered(
        return_df: bool = False,
    ) -> Generator[tuple] | pd.DataFrame:
    """
    Filtered textmining dataset of the DISEASES database from Jensen Lab.

    Warning: The "textmining" datasets are enormous!

    Args:
        return_df:
            Return a pandas data frame.
    """

    return diseases_general(
        data_origin = 'textmining',
        filtered = True,
        return_df = return_df,
    )


def knowledge_filtered(
        return_df: bool = False,
    ) -> Generator[tuple] | pd.DataFrame:
    """
    Filtered knowledge dataset of the DISEASES database from Jensen Lab.

    Args:
        return_df:
            Return a pandas data frame.
    """

    return diseases_general(
        data_origin = 'knowledge',
        filtered = True,
        return_df = return_df,
    )


def knowledge_full(
        return_df: bool = False,
    ) -> Generator[tuple] | pd.DataFrame:
    """
    Full knowledge dataset of the DISEASES database from Jensen Lab.

    Args:
        return_df:
            Return a pandas data frame.
    """

    return diseases_general(
        data_origin = 'knowledge',
        filtered = False,
        return_df = return_df,
    )


def experiments_filtered(
        return_df: bool = False,
    ) -> Generator[tuple] | pd.DataFrame:
    """
    Filtered experiments dataset of the DISEASES database from Jensen Lab.

    Args:
        return_df:
            Return a pandas data frame.
    """

    return diseases_general(
        data_origin = 'experiments',
        filtered = True,
        return_df = return_df,
    )


def experiments_full(
        return_df: bool = False,
    ) -> Generator[tuple] | pd.DataFrame:
    """
    Full experiments dataset of the DISEASES database from Jensen Lab.

    Args:
        return_df:
            Return a pandas data frame.
    """

    return diseases_general(
        data_origin = 'experiments',
        filtered = False,
        return_df = return_df,
    )
