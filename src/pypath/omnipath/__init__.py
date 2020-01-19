#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module.
#  Provides a high level interface for managing builds of the
#  OmniPath databases.
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


import os
import copy

from pypath import session_mod as _session_mod
from pypath import settings as _settings_mod
from pypath.omnipath import app as _app_mod


_logger = _session_mod.Logger(name = 'omnipath.init')
_log = _logger._log


_log('Welcome to the OmniPath database manager app.')


def init(**kwargs):

    param = (
        copy.deepcopy(globals()['OP_DB_ARGS'])
            if 'OP_DB_ARGS' in globals() else
        {}
    )

    _log(
        'You can customize the database building process by '
        'setting parameters in the `OP_DB_ARGS` global variable '
        'or by calling `init` again with keyword arguments or after '
        'setting values in the `pypath.settings` module.'
    )

    param.update(kwargs)

    globals()['data'] = _app_mod.Database(**param)


init()
