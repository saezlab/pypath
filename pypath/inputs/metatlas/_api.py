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

"""
Client for Metabolic Atlas REST API.

Metabolic Atlas provides access to Genome-Scale Metabolic Models (GEMs)
from various organisms and tissues.
"""

from __future__ import annotations

import json

import pypath.share.curl as curl
import pypath.resources.urls as urls

from ._common import _log, REQ_HEADERS, _df_or_records
from ._records import MetatlasModel

__all__ = [
    'metatlas_models',
    'metatlas_integrated_models',
]


def metatlas_models(dataframe: bool = False):
    """
    Retrieves all repository models from Metabolic Atlas.

    Metabolic Atlas contains 360+ metabolic models from various
    organisms and tissues.

    Args:
        dataframe: Return a pandas DataFrame instead of a list.

    Returns:
        List of MetatlasModel named tuples, or DataFrame if requested.
    """

    _log('Downloading repository models from Metabolic Atlas.')

    url = urls.urls['metatlas']['models']
    c = curl.Curl(url, silent=False, large=False, req_headers=REQ_HEADERS)

    if c.result is None:
        _log('Failed to download Metabolic Atlas models.')
        return [] if not dataframe else None

    data = json.loads(c.result)
    records = [_parse_model(model) for model in data]

    return _df_or_records(records, dataframe)


def metatlas_integrated_models(dataframe: bool = False):
    """
    Retrieves integrated models from Metabolic Atlas.

    Integrated models are 7 curated models with enhanced API access
    and additional features like compartment maps and subsystem data.

    Args:
        dataframe: Return a pandas DataFrame instead of a list.

    Returns:
        List of MetatlasModel named tuples, or DataFrame if requested.
    """

    _log('Downloading integrated models from Metabolic Atlas.')

    url = urls.urls['metatlas']['integrated']
    c = curl.Curl(url, silent=False, large=False, req_headers=REQ_HEADERS)

    if c.result is None:
        _log('Failed to download Metabolic Atlas integrated models.')
        return [] if not dataframe else None

    data = json.loads(c.result)
    records = [_parse_model(model) for model in data]

    return _df_or_records(records, dataframe)


def _parse_model(model: dict) -> MetatlasModel:
    """
    Parses a model dictionary from the API into a MetatlasModel.

    Handles both repository models and integrated models which have
    different response structures.

    Args:
        model: Dictionary from the API response.

    Returns:
        MetatlasModel named tuple.
    """

    sample = model.get('sample', {}) or {}
    gemodelset = model.get('gemodelset', {}) or {}

    return MetatlasModel(
        id=model.get('id', ''),
        name=model.get('full_name', '') or gemodelset.get('name', ''),
        short_name=model.get('short_name', ''),
        organism=sample.get('organism', ''),
        tissue=sample.get('tissue', ''),
        cell_type=sample.get('cell_type', ''),
        condition=model.get('condition', ''),
        year=model.get('year') or model.get('date', ''),
        reaction_count=model.get('reaction_count'),
        metabolite_count=model.get('metabolite_count'),
        gene_count=model.get('gene_count'),
        files=tuple(model.get('files', [])) if model.get('files') else (),
    )
