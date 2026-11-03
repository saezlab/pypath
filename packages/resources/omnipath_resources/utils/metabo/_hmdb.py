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

"""
Utilities for HMDB identifier normalisation.
"""

__all__ = ['normalise_hmdb']

_HMDB_DIGITS = 7  # current zero-padded width of the numeric part


def normalise_hmdb(hmdb_id: str) -> str:
    """
    Normalise an HMDB identifier to the current 7-digit zero-padded format.

    HMDB IDs were historically distributed in a 5-digit format
    (e.g. ``HMDB00001``).  The current format uses 7 digits
    (e.g. ``HMDB0000001``).  IDs already in the new format are returned
    unchanged.  Non-HMDB strings are returned as-is.

    Args:
        hmdb_id: HMDB identifier string (e.g. ``'HMDB00001'`` or
            ``'HMDB0000001'``).

    Returns:
        HMDB ID with 7-digit zero-padded numeric part
        (e.g. ``'HMDB0000001'``), or the original string unchanged if
        it does not start with ``'HMDB'``.
    """

    if not hmdb_id.upper().startswith('HMDB'):
        return hmdb_id

    return 'HMDB' + hmdb_id[4:].zfill(_HMDB_DIGITS)
