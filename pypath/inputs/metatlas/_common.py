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
Common utilities and logger for Metabolic Atlas input module.
"""

from __future__ import annotations

import pypath.share.session as session

__all__ = ['_log', 'REQ_HEADERS']

_log = session.Logger(name='metatlas_input')._log

REQ_HEADERS = [
    'Accept: application/json',
    'User-Agent: pypath (https://pypath.omnipathdb.org)',
]
