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

import pypath.data as _data


def synonyms(resource: str) -> dict[str, str]:
    """
    ID type synonyms for one resource.

    Args:
        resource:
            Name of the resource.

    Returns:
        A dict with synonyms as keys and resource ID types as values.
    """

    syn = _data.common_load(f'idtypes/{resource.lower()}')

    if isinstance(syn, list):

        syn = {x: x for x in syn}

    return syn
