#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2018
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://www.ebi.ac.uk/~denes
#

# XXX: Missing module docstring

# try:
#import __main__ as ext
from __future__ import print_function

import pypath.main as main
import pypath._version as _version
#import descriptions
#import pypath.common as common

__version__ = _version.__version__
PyPath = main.PyPath # This is an alias
