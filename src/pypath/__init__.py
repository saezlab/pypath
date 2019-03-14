#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2019
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#                  Nicolàs Palacio
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

"""
``pypath`` is a module primarily for building molecular interaction networks
but also with several submodules for accessing, preprocessing and serving
data from various resources.
"""

import pypath._version as _version_mod
import pypath.session as _session_mod

__version__ = _version.__version__
__author__ = _version.__author__

_session_mod.new_session()
session = _session_mod.get_session()
