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

from typing import Optional

import os

import pypath.share.settings as settings
import pypath.share.session as session

_logger = session.Logger(name = 'credentials')
_log = _logger._log


def credentials(
        *args: tuple[str, str],
        resource: Optional[str] = None,
        from_file: Optional[str] = None,
        **kwargs: dict[str, str],
    ) -> dict:
    """
    Credentials required for restricted access resources.

    Args
        args:
            Two strings: a user name and password. If only one provided, it
            is assumed to be a user name; if more provided, apart from the
            first two, the rest will be ignored.
        resource:
            Name of the resource. If the key `<resource>_credentials`
            exists in the module settings, its value will be returned as
            credentials.
        from_file:
            Path to a file or name of a file that is located in the module's
            default secrets directory.
        kwargs:
            Custom key-value pairs, will be returned unchanged. This is the
            way to explicitely provide user and password, and any further
            fields.

    Returns
        A dictionary with the credentials. Raises RuntimeError if credentials
        not provided by any of the available ways.
    """

    fields = ('user', 'passwd')
    kwargs.update(dict(zip(fields, args)))
    kwargs = dict(it for it in kwargs.items() if it[1] is not None)

    if all(f in kwargs for f in fields):

        credentials = kwargs

    else:

        settings_key = f'{resource.lower()}_credentials'
        credentials = settings.get(settings_key)

        if not credentials:

            secrets_fname = from_file or settings_key

            if not os.path.exists(secrets_fname):

                secrets_fname = os.path.join(
                    settings.get('secrets_dir'),
                    secrets_fname,
                )

            if os.path.exists(secrets_fname):

                _log(
                    f'Reading credentials for `{resource}` '
                    f'from file `{secrets_fname}`.'
                )

                with open(secrets_fname, 'r') as fp:

                    lines = fp.read().strip().split(os.linesep)

                keys, values = tuple(zip(*(
                    ([None] + l.split(':', maxsplit = 1))[-2:]
                    for l in lines
                )))

                keys = keys if all(keys) else fields
                credentials = dict(zip(keys, values))
                credentials.update(kwargs)

        else:

            _log(f'`{resource}` credentials provided by `settings`.')

        if not credentials:

            msg = f'Failed to obtain credentials for resource `{resource}`'
            _log(msg)

            raise RuntimeError(msg)

    return credentials
