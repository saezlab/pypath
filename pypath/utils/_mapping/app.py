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


def init(**kwargs):
    """
    Create a new `Mapper` instance under the `mapper` attribute of this
    module.

    Returns
        None.
    """

    if 'mapper' in globals():

        globals()['mapper'].__del__()

    globals()['mapper'] = Mapper(**kwargs)


def get_mapper(**kwargs):
    """
    The module under its `mapper` attribute has an instance of the `Mapper`
    object, which manages the ID translations. This function creates the
    instance if does not exist and returns it.

    Returns
        A Mapper object.
    """

    if 'mapper' not in globals():

        init(**kwargs)

    return globals()['mapper']
