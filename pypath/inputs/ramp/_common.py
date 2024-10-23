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

import pprint

import pypath.share.session as session

__all__ = ['_log', '_show_tables']

_log = session.Logger(name = 'ramp_input')._log


def _show_tables(tables: dict[str, list[str]]) -> None:
    """
    Show the tables of the RaMP database from SQL dump.
    """

    pprint.pprint(tables)


