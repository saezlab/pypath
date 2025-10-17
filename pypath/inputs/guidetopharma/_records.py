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

from typing import NamedTuple

__all__ = [
    'G2PInteraction',
    'G2PLigand',
    'G2PTarget',
    'clean_dict',
]

def clean_dict(d: dict) -> dict:
    """Return a new dict with empty strings converted to None."""
    return {k: (None if v == '' else v) for k, v in d.items()}


class G2PLigand(NamedTuple):
    """Guide to Pharmacology ligand record."""
    ligand_id: str | None = None
    name: str | None = None
    uniprot: str | None = None
    pubchem: str | None = None
    iupac: str | None = None
    chembl: str | None = None
    smiles: str | None = None
    inchi: str | None = None
    organism: str | None = None
    entity_type: str | None = 'ligand'
    subtype: str | None = None
    role: str | None = None


class G2PTarget(NamedTuple):
    """Guide to Pharmacology target record."""
    target_id: str | None = None
    uniprot: str | None = None
    symbol: str | None = None
    entrez: str | None = None
    ensembl: str | None = None
    refseq: str | None = None
    refseqp: str | None = None
    family: str | None = None
    target_type: str | None = None
    organism: str | None = None
    entity_type: str | None = 'target'
    role: str | None = None


class G2PInteraction(NamedTuple):
    """Guide to Pharmacology interaction record."""
    ligand: G2PLigand | None
    target: G2PTarget
    action: str | None = None
    action_type: str | None = None
    is_stimulation: bool | None = None
    is_inhibition: bool | None = None
    endogenous: bool | None = None
    affinity_high: str | None = None
    affinity_low: str | None = None
    affinity_median: str | None = None
    affinity_units: str | None = None
    primary_target: str | None = None
    pubmed: str | None = None
