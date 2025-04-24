#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright 2014-2025
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

from typing import Literal
from collections.abc import Generator

import csv

from pypath.share import curl
from pypath.resources import urls

__all__ = [
    'table',
]


TABLES = Literal[
    'ligands',
    'interactions',
    'targets_and_families',
    'ligand_id_mapping',
    'ligand_physchem_properties',
    'endogenous_ligand_pairings_all',
    'endogenous_ligand_detailed',
    'approved_drug_primary_target_interactions',
    'approved_drug_detailed_interactions',
    'peptides',
    'GtP_to_HGNC_mapping',
    'GtP_to_UniProt_mapping',
]


def table(name: TABLES) -> Generator[dict]:
    """
    Downloads a table from Guide2Pharma.

    Args:
        name
            The name of the table to download.

    Returns:
        Generator yielding the table as named tuples.
    """

    url = urls.urls['gtp']['url'] % name

    c = curl.Curl(url, silent = False, large = True, encoding = 'utf-8')

    g2p_version = next(c.result).strip()

    yield from csv.DictReader(c.result)
