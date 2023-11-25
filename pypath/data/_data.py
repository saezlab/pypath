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

import itertools
import functools

from pypath_common import data as _data

_METHODS = ('load', 'path', 'builtins')
_MODULES = {
    'pypath': '',
    'pypath_common': 'common_',
    'pypath.resources': 'resources_',
}

for method, (module, prefix) in itertools.product(_METHODS, _MODULES.items()):

    globals()[prefix + method] = functools.partial(
        getattr(_data, method),
        module = module,
    )
