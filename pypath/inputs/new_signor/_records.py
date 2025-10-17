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

from typing import NamedTuple

__all__ = [
    'SignorComplex',
    'SignorProteinFamily',
    'SignorPhenotype',
    'SignorStimulus',
]


class SignorComplex(NamedTuple):
    """SIGNOR complex record."""
    complex_id: str | None = None
    name: str | None = None
    components: list[str] | None = None


class SignorProteinFamily(NamedTuple):
    """SIGNOR protein family record."""
    family_id: str | None = None
    name: str | None = None
    members: list[str] | None = None


class SignorPhenotype(NamedTuple):
    """SIGNOR phenotype record."""
    phenotype_id: str | None = None
    name: str | None = None
    description: str | None = None


class SignorStimulus(NamedTuple):
    """SIGNOR stimulus record."""
    stimulus_id: str | None = None
    name: str | None = None
    description: str | None = None
