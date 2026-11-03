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
Raw download functions for iMM1415.

The primary source is the BiGG JSON model downloaded from BiGG Models.
Unlike Recon3D, iMM1415 has no supplementary MATLAB file from VMH.
"""

from __future__ import annotations

__all__ = ['imm1415_raw']

import json

import pypath.resources.urls as urls
import pypath.share.curl as curl

from ._common import _log


def imm1415_raw() -> dict:
    """
    Download and parse the iMM1415 model from BiGG in JSON format.

    The JSON contains ``metabolites``, ``reactions``, ``genes``, and
    ``compartments`` keys.  Metabolites and genes carry cross-reference
    annotations (HMDB, ChEBI, KEGG, Entrez, etc.) in their ``annotation``
    fields.

    Downloaded once per session; pypath's curl cache prevents repeated
    network requests.

    Returns:
        Parsed JSON as a Python dict.
    """

    url = urls.urls['imm1415']['bigg_json']
    _log(f'Downloading iMM1415 JSON from BiGG: {url}')
    c = curl.Curl(url, silent=False, large=True)

    if c.result is None:
        _log('Failed to download iMM1415 JSON from BiGG.')
        return {}

    if isinstance(c.result, str):
        return json.loads(c.result)

    return json.loads(''.join(c.result))
