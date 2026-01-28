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

from collections.abc import Generator
import json

import pypath.share.curl as curl
import pypath.resources.urls as urls

from ._common import _log, REQ_HEADERS
from ._records import MetatlasModel

__all__ = [
    'metatlas_models',
    'metatlas_integrated_models',
    'metatlas_model_files',
]


def metatlas_models() -> Generator[MetatlasModel, None, None]:
    """
    Retrieves all repository models from Metabolic Atlas.

    Metabolic Atlas contains 360+ metabolic models from various
    organisms and tissues.

    Yields:
        MetatlasModel named tuples with model metadata.
    """

    _log('Downloading repository models from Metabolic Atlas.')

    url = urls.urls['metatlas']['models']
    c = curl.Curl(url, silent=False, large=False, req_headers=REQ_HEADERS)

    if c.result is None:
        _log('Failed to download Metabolic Atlas models.')
        return

    data = json.loads(c.result)

    for model in data:
        yield _parse_model(model)


def metatlas_integrated_models() -> Generator[MetatlasModel, None, None]:
    """
    Retrieves integrated models from Metabolic Atlas.

    Integrated models are 7 curated models with enhanced API access
    and additional features like compartment maps and subsystem data.

    Yields:
        MetatlasModel named tuples with model metadata.
    """

    _log('Downloading integrated models from Metabolic Atlas.')

    url = urls.urls['metatlas']['integrated']
    c = curl.Curl(url, silent=False, large=False, req_headers=REQ_HEADERS)

    if c.result is None:
        _log('Failed to download Metabolic Atlas integrated models.')
        return

    data = json.loads(c.result)

    for model in data:
        yield _parse_model(model)


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

    # Integrated models use full_name, repository models use gemodelset.name
    name = model.get('full_name', '') or gemodelset.get('name', '')

    # Integrated models use date, repository models use year
    year = model.get('year') or model.get('date', '')

    return MetatlasModel(
        id=model.get('id', ''),
        name=name,
        short_name=model.get('short_name', ''),
        organism=sample.get('organism', ''),
        tissue=sample.get('tissue', ''),
        cell_type=sample.get('cell_type', ''),
        condition=model.get('condition', ''),
        year=year,
        reaction_count=model.get('reaction_count'),
        metabolite_count=model.get('metabolite_count'),
        gene_count=model.get('gene_count'),
        files=tuple(model.get('files', [])) if model.get('files') else (),
    )


def metatlas_model_files(model_id: str) -> list[dict]:
    """
    Retrieves file list for a specific model.

    Args:
        model_id: The model ID (e.g., 'Human-GEM').

    Returns:
        List of file dictionaries with path, name, and type info.
    """

    _log(f'Getting file list for model {model_id}.')

    url = urls.urls['metatlas']['file'] % model_id
    c = curl.Curl(url, silent=False, large=False, req_headers=REQ_HEADERS)

    if c.result is None:
        _log(f'Failed to get files for model {model_id}.')
        return []

    data = json.loads(c.result)
    return data.get('files', [])


def _download_model_file(path: str) -> str | None:
    """
    Downloads a file from a model repository.

    Args:
        path: The file path from the model files list.

    Returns:
        File content as string, or None if download failed.
    """

    url = urls.urls['metatlas']['file'] % path
    c = curl.Curl(url, silent=False, large=True, req_headers=REQ_HEADERS)

    if c.result is None:
        _log(f'Failed to download file {path}.')
        return None

    return c.result
