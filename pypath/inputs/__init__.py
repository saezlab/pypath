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

import importlib


import pypath.share.session as session

_logger = session.Logger(name = 'inputs')
_log = _logger._log


def get_method(module_name, method_name = None):
    """
    Retrieves a method from a submodule of this module (``inputs``) by its
    name. E.g. for ``'uniprot_db.all_uniprots'`` it returns the ``all_uniprots``
    method from the ``pypath.inputs.uniprot_db`` module.
    """

    _log('Selecting input method (step 1): module `%s`, method `%s`.' % (
            module_name,
            method_name,
        )
    )

    if callable(module_name):

        return module_name

    if not method_name:

        module_method = module_name.rsplit('.', maxsplit = 1)
        method_name = module_method[-1]
        module_name = module_method[-2] if len(module_method) > 1 else 'main'

    module_name = module_name.rsplit('.', maxsplit = 1)[-1]
    module_name = 'pypath.inputs.%s' % module_name

    _log('Selecting input method (step 2): module `%s`, method `%s`.' % (
            module_name,
            method_name,
        )
    )

    try:

        _log('Importing module `%s`.' % module_name)
        mod = importlib.import_module(module_name)

    except:

        session.log().msg(
            msg = 'Could not import module `%s`.' % module_name,
            label = 'inputs',
        )

    try:

        method = getattr(mod, method_name)

        return method

    except:

        session.log().msg(
            msg = 'Could not find method `%s` in module `%s`.' % (
                method_name,
                module_name,
            ),
            label = 'inputs',
        )
