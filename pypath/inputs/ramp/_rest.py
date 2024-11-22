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

import json

import pypath.resources.urls as urls
import pypath.share.curl as curl

__all__ = ['ramp_id_types']


def ramp_id_types(
        entity_type: Literal['gene', 'compound'] | None = None,
    ) -> set[str]:
    """
    List the identifier types of the RaMP database.
    """

    entity_types = {
        'compound': 'Metabolites',
        'gene': 'Genes/Proteins',
    }

    url = urls.urls['ramp']['api'] % 'id-types'
    c = curl.Curl(url, silent = True, large = False)

    return {
        id_type.strip()
        for i in json.loads(c.result)['data']
        if not entity_type or i['analyteType'] == entity_types[entity_type]
        for id_type in i['idTypes'].split(',')
    }
