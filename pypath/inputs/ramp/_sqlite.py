#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright 2014-2024
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

import os
import sqlite3

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.cache as cache
from ._common import _log


def ramp_sqlite(version: str = '2.5.4') -> sqlite3.Connection:
    """
    Download the RaMP database in SQLite format.

    SQLite format builds are avialable in RaMP's GitHub repo:
        https://github.com/ncats/RaMP-DB/tree/main/db
    This function downloads and opens one of them.

    Args
        version:
            The version of the RaMP database to download.

    Returns
        A SQLite database connection.
    """

    sqlite_path = os.path.join(
        cache.get_cachedir(),
        f'RaMP_SQLite_{version}.sqlite',
    )
    exists = os.path.exists(sqlite_path)
    _log(f'RaMP SQLite path: `{sqlite_path}`; exists: {exists}.')

    if not exists:

        url = urls.urls['ramp']['github_sqlite'] % version
        c = curl.Curl(url, large = True, silent = False, compr = 'gz')

        with open(sqlite_path, 'wb') as fp:

            while chunk := c.gzfile.read(1024 ** 2):

                fp.write(chunk)

    if os.path.exists(sqlite_path):

        _log(f'Opening SQLite: `{sqlite_path}`')
        con = sqlite3.connect(sqlite_path)

        return con
