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
Visualisation of metabolite structures.
"""

from typing import Iterable

import pypath.share.common as common
import pypath.utils.mapping as mapping
import pypath.inputs.hmdb.visual as hmdb_visual


def show(metabolite: str | Iterable, id_type: str | None = None) -> None:
    """
    Show the 2D structure(s) of one or more metabolite(s).

    Args:
        metabolite:
            One or more metabolite identifier(s).
        id_type:
            Identifier type of the metabolite(s). If not provided, HMDB
            identifier is assumed.
    """

    if isinstance(metabolite, common.LIST_LIKE):

        for m in metabolite:

            show(m, id_type = id_type)

        return

    metabolite = (
        (metabolite,)
            if metabolite.startswith('HMDB') else
        mapping.map_name(metabolite, id_type, 'hmdb')
    )

    for m in metabolite:

        hmdb_visual.show_structure(m)
