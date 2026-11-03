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

from typing import Generator

import collections

import pypath.share.curl as curl
import pypath.resources.urls as urls
from pypath.inputs.offsides import _sides_base


def twosides_interactions() -> Generator[tuple]:
    """
    Drug-drug interaction data from the TwoSIDES database.

    Retrieves drug-drug interaction safety signals mined
    from the FDA's Adverse Event Reporting System.

    Yields:
        Tuples of drug-drug interactions.
    """

    return _sides_base(
        url_key = 'twosides',
        fields = (
            'drug1_rxnorm',
            'drug1',
            'drug2_rxnorm',
            'drug2',
            'condition_meddra',
            'condition',
            'prr',
            'pee_error',
            'mean_reporting_frequency',
        ),
        indices = (0, 1, 2, 3, 4, 5, 10, 11, 12),
        record_name = 'TwosidesInteraction',
    )
