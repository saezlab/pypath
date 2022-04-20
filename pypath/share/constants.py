#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2022
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: Dénes Türei (turei.denes@gmail.com)
#           Nicolàs Palacio
#           Sebastian Lobentanzer
#           Erva Ulusoy
#           Olga Ivanova
#           Ahmet Rifaioglu
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

NO_VALUE = 'PYPATH_NO_VALUE'
GLOM_ERROR = 'PYPATH_GLOM_ERROR'
CURSOR_UP_ONE = '\x1b[1A'
ERASE_LINE = '\x1b[2K'
NOT_ORGANISM_SPECIFIC = -1
BOOLEAN_TRUE = frozenset(('1', 'yes', 'true'))
BOOLEAN_FALSE = frozenset(('0', 'no','false'))
BOOLEAN_VALUES = BOOLEAN_TRUE.union(BOOLEAN_FALSE)
