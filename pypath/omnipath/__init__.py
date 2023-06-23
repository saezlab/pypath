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

import os
import copy

from pypath.share import session as _session_mod
from pypath.share import settings as _settings_mod
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
        'setting values in the `pypath.share.settings` module.'
        'You can define your own databases following the examples in '
        '`pypath/omnipath/databases/builtins.json`. '
        'See more details in `pypath/omnipath/databases/db_template.json` '
        'and the code of method '
        '`pypath.omnipath.app.DatabaseManager.load_dataset`.'
    )

    param.update(kwargs)

    globals()['db'] = _app_mod.DatabaseManager(**param)


init()
