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

import tarfile
import sqlite3
import shutil
import os

import pypath.share.cache as cache
from pypath.share import curl, cache
from pypath.resources import urls
from ._common import _log, _show_tables


def _chembl_sqlite_path(version: int = 36) -> str:

    return os.path.join(
        cache.get_cachedir(),
        f'ChEMBL_SQLite_{version:02}.sqlite',
    )

def chembl_sqlite(
        version: int = 36,
        connect: bool = True,
    ) -> sqlite3.Connection | str:

    """
    Download the ChEMBL database in SQLite format.

    SQLite format builds are avialable in ChEMBL's FTP site:
        https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases
    This function downloads and opens one of them.

    Args
        version:
            The version of the ChEMBL database to download.

    Returns
        A SQLite database connection.
    """

    sqlite_path = _chembl_sqlite_path(version)
    exists = os.path.exists(sqlite_path)
    _log(f'ChEMBL SQLite path: `{sqlite_path}`; exists: {exists}.')

    if not exists:

        path_in_tar = (
            'chembl_%02i/chembl_%02i_sqlite/chembl_%02i.db' %
            ((version, ) * 3)
        )
        url = urls.urls['chembl']['sqlite'] % (version, version)
        c = curl.Curl(
            url,
            large = True,
            silent = False,
            files_needed = [path_in_tar],
            slow = True
        )
        with open(sqlite_path, 'wb') as fp:

            while chunk := c.result[path_in_tar].read(1024 ** 2):

                fp.write(chunk)

    if connect:

        _log(f'Opening SQLite: `{sqlite_path}`')
        con = sqlite3.connect(sqlite_path)

        return con

    return sqlite_path
