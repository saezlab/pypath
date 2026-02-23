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
Raw download functions for Recon3D.

Two sources are available:

- **BiGG** (primary): JSON model from BiGG Models database.  Contains
  stoichiometry, gene-reaction rules, and metabolite/gene cross-references
  including HMDB and ChEBI annotations.
- **VMH** (supplementary): MATLAB ``.mat`` file from Virtual Metabolic Human.
  Provides an independent set of HMDB annotations for metabolites, accessed
  via :func:`gem_matlab_extract` (the Python equivalent of OmnipathR's
  ``gem_matlab_tibble``).

The MATLAB parsing requires ``scipy``.
"""

from __future__ import annotations

__all__ = [
    'recon3d_raw',
    'recon3d_raw_vmh',
    'gem_matlab_extract',
]

import json

import pypath.resources.urls as urls
import pypath.share.curl as curl

from ._common import _log


def recon3d_raw() -> dict:
    """
    Download and parse the Recon3D model from BiGG in JSON format.

    The JSON contains ``metabolites``, ``reactions``, ``genes``, and
    ``compartments`` keys.  Metabolites and genes carry cross-reference
    annotations (HMDB, ChEBI, KEGG, Entrez, etc.) in their ``annotation``
    fields.

    Downloaded once per session; pypath's curl cache prevents repeated
    network requests.

    Returns:
        Parsed JSON as a Python dict.
    """

    url = urls.urls['recon3d']['bigg_json']
    _log(f'Downloading Recon3D JSON from BiGG: {url}')
    c = curl.Curl(url, silent=False, large=True)

    if c.result is None:
        _log('Failed to download Recon3D JSON from BiGG.')
        return {}

    # With large=True, c.result is a file-like object; use json.load not json.loads.
    try:
        return json.load(c.result)
    except (TypeError, AttributeError):
        # Fallback: result may already be a string in some curl versions.
        return json.loads(c.result)


def gem_matlab_extract(mat: dict, *fields: str) -> dict[str, list]:
    """
    Extract named fields from a MATLAB model loaded by ``scipy.io.loadmat``.

    Python equivalent of OmnipathR's ``gem_matlab_tibble()``.  The MATLAB
    object loaded by scipy is a deeply nested numpy array structure; this
    function navigates the nesting to return plain Python lists.

    The top-level model key (e.g. ``'Recon3D'``, ``'ihuman'``, ``'model'``)
    is auto-detected as the first key that does not start with ``'_'``.

    Args:
        mat:
            Dict returned by ``scipy.io.loadmat()``.
        *fields:
            MATLAB field names to extract (e.g. ``'mets'``, ``'metHMDBID'``).

    Returns:
        Dict mapping each field name to a list of extracted values.  Scalar
        entries become Python scalars; vector entries become lists of strings.
    """

    # Auto-detect the top-level model key (skip scipy metadata keys).
    model_key = next((k for k in mat if not k.startswith('_')), None)

    if model_key is None:
        _log('gem_matlab_extract: no model key found in MATLAB object.')
        return {f: [] for f in fields}

    model = mat[model_key][0, 0]
    result: dict[str, list] = {}

    for field in fields:

        if field not in model.dtype.names:
            _log(f'gem_matlab_extract: field {field!r} not found in MATLAB model.')
            result[field] = []
            continue

        raw = model[field]  # typically shape (n, 1)
        values: list = []

        for i in range(raw.shape[0]):

            elem = raw[i, 0]

            # Scalar: all dimensions are 1 (matches R's `dim %>% equals(1L) %>% all`).
            if elem.size == 1:
                val = elem.flat[0]
                values.append(
                    val.decode('utf-8', errors='replace')
                    if isinstance(val, bytes)
                    else val
                )

            else:
                # Vector/array: flatten to a list of strings.
                values.append([
                    v.decode('utf-8', errors='replace') if isinstance(v, bytes) else v
                    for v in elem.flat
                ])

        result[field] = values

    return result


def recon3d_raw_vmh() -> dict:
    """
    Download Recon3D from Virtual Metabolic Human (VMH) as a MATLAB ``.mat``
    file and parse it with ``scipy.io.loadmat``.

    Used to obtain supplementary metabolite annotations (e.g. HMDB IDs) that
    are not available or less complete in the BiGG JSON.  Requires ``scipy``.

    Returns:
        Dict returned by ``scipy.io.loadmat``.  Pass to
        :func:`gem_matlab_extract` to extract specific fields.
    """

    try:
        import scipy.io

    except ImportError as exc:
        raise ImportError(
            'scipy is required for MATLAB parsing. '
            'Install it with: pip install scipy'
        ) from exc

    url = urls.urls['recon3d']['vmh_mat']
    _log(f'Downloading Recon3D MATLAB from VMH: {url}')
    c = curl.Curl(url, silent=False, large=True)

    if c.result is None:
        _log('Failed to download Recon3D MATLAB from VMH.')
        return {}

    return scipy.io.loadmat(c.fileobj.name)
