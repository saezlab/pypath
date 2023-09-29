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


def get_controller(reload = False):
    """
    Returns the resource controller.

    The controller instantiated only once by default, at the first access
    attempt. The instance is stored in the module and provided on demand.
    """

    from . import controller as _controller_mod

    if '_controller' not in globals() or reload:

        globals()['_controller'] = _controller_mod.ResourceController()

    return globals()['_controller']
