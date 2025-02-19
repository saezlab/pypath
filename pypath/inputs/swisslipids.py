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

"""
Access SwissLipids datasets.

https://swisslipids.org/#/downloads
"""

from typing import Literal
from collections.abc import Generator

import re
import csv
import functools

import pandas as pd

from pypath_common import _misc
from pypath.share import curl as _curl
from pypath.resources import urls as _urls
from pypath.share import session as _session

_logger = _session.Logger(name = 'swisslipids_input')
_log = _logger._log


DATASETS = Literal[
    'enzymes',
    'evidences',
    'lipids',
    'go',
    'lipids2uniprot',
    'tissues'
]
__all__ = [f'swisslipids_{_}' for _ in DATASETS.__args__]


def _swisslipids(
        dataset: DATASETS,
        return_df: bool = False,
    ) -> Generator[dict[str, str]] | pd.DataFrame:
    """
    Access a SwissLipids dataset.

    Args:
        dataset:
            Name of the dataset, as shown at
            https://swisslipids.org/#/downloads. Available datasets are:
            - enzymes
            - evidences
            - lipids
            - go
            - lipids2uniprot
            - tissues
        return_df:
            Return a pandas data frame. The "lipids" and "lipids2uniprot"
            datasets are large and require up to 1.1 GB of memory if loaded
            into a data frame.

    Returns:
        A generator yielding dicts, each representing a record in the dataset;
        or a pandas data frame if `return_df` is True.
    """

    _log(f'Loading SwissLipids dataset `{dataset}`.')
    url = _urls.urls['swisslipids']['url'] % dataset

    c = _curl.Curl(
        url,
        large = True,
        silent = False,
        encoding = 'latin-1',
        compr = None if dataset == 'go' else 'gz',
    )

    lines = (re.sub('\t\\s+', '\t', ll).strip() for ll in c.result)

    if return_df:

        _log(f'Loading SwissLipids `{dataset}` dataset into data frame.')
        lines = (l.split('\t') for l in lines)
        cols = next(lines)
        df = pd.DataFrame(lines, columns=cols) #c._gzfile_mode_r
        _log(
            f'SwissLipids `{dataset}` data frame '
            f'memory usage: `{_misc.df_memory_usage(df)}`.'
        )

        return df

    else:

        _log(f'Iterating records in SwissLipids `{dataset}` dataset.')
        return csv.DictReader(lines, delimiter = '\t')



for _dataset in DATASETS.__args__:

    _func_name = f'swisslipids_{_dataset}'
    globals()[_func_name] = functools.partial(_swisslipids, dataset = _dataset)
