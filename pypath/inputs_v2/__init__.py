#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright 2014-2025
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

"""
Input modules v2 - aligned with the new internal data representation.

These modules emit Entity records using the schema defined in
omnipath_build.utils.silver_schema, replacing the legacy input modules
in pypath.inputs.
"""

import importlib


import pypath.share.session as session

_logger = session.Logger(name='inputs_v2')
_log = _logger._log


def get_method(module_name, method_name=None):
    """
    Retrieves a method from a submodule of this module (``inputs_v2``) by its
    name. E.g. for ``'uniprot.uniprot_entities'`` it returns the
    ``uniprot_entities`` method from the ``pypath.inputs_v2.uniprot`` module.

    This function works similarly to the v1 inputs.get_method but is designed
    for the new v2 modules that emit Entity records according to the new schema.
    """

    _log(
        'Selecting input method v2 (step 1): module `%s`, method `%s`.'
        % (module_name, method_name)
    )

    if callable(module_name):
        return module_name

    if not method_name:
        module_method = module_name.rsplit('.', maxsplit=1)
        method_name = module_method[-1]
        module_name = module_method[-2] if len(module_method) > 1 else 'main'

    module_name = module_name.rsplit('.', maxsplit=1)[-1]
    module_name = 'pypath.inputs_v2.%s' % module_name

    _log(
        'Selecting input method v2 (step 2): module `%s`, method `%s`.'
        % (module_name, method_name)
    )

    try:
        _log('Importing module `%s`.' % module_name)
        mod = importlib.import_module(module_name)

    except Exception as e:
        _logger._logger.msg(
            msg='Could not import module `%s`: %s' % (module_name, str(e)),
            label='inputs_v2',
        )
        raise

    try:
        method = getattr(mod, method_name)
        return method

    except AttributeError as e:
        _logger._logger.msg(
            msg='Could not find method `%s` in module `%s`: %s'
            % (method_name, module_name, str(e)),
            label='inputs_v2',
        )
        raise
