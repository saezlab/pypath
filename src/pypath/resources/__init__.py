#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2020
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#                  Nicolàs Palacio
#                  Olga Ivanova
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

from pypath.resources import controller as _controller_mod


def get_controller(reload = False):
    """
    Returns the resource controller.

    The controller instantiated only once by default, at the first access
    attempt. The instance is stored in the module and provided on demand.
    """

    if '_controller' not in globals() or reload:

        globals()['_controller'] = _controller_mod.ResourceController()

    return globals()['_controller']
