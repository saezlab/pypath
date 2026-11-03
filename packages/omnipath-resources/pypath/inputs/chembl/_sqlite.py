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

import sqlite3
from typing import Callable
from pathlib import Path

import pypath.share.cache as cache
from pypath.share import curl, cache
from pypath.formats import sqlite as _sqlite
from pypath.resources import urls

from ._common import _log, _show_tables

def chembl_sqlite(
        version: int = 36,
        connect: bool = True,
    ) -> sqlite3.Connection | Path:

    """
    Download the ChEMBL database in SQLite format.

    SQLite format builds are available at ChEMBL's FTP site:
    https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases
    This function downloads, caches, and opens one of them.

    Args
        version:
            The version of the ChEMBL database to download.

    Returns:
        A SQLite database connection to the ChEMBL database, or the path to
        the database file.
    """

    path_in_tar = (
        f'chembl_{version:02d}/chembl_{version:02d}_sqlite/chembl_{version:02d}.db'
        )
    url = urls.urls['chembl']['sqlite'] % (version, version)

    def _chembl_download() -> curl.Curl:
        """Callback to download the tarballed ChEMBL database."""
        return curl.Curl(
            url,
            large = True,
            silent = False,
            files_needed = [path_in_tar],
            slow = True,
        )

    def _chembl_extractor(curl_obj: curl.Curl) -> str:
        """Extractor to get the SQLite file from the tarball."""
        return curl_obj.result[path_in_tar]

    return _sqlite.download_sqlite(
        download_callback= _chembl_download,
        extractor=_chembl_extractor,
        database = 'ChEMBL',
        version = f'{version:02}',
        connect = connect,
    )

def chembl_sqlite_raw() -> dict[str, list[tuple]]:
    """Extracts all tables and data from the ChEMBL SQLite database.

    Returns:
        A dictionary where keys are table names and values are lists of
        tuples representing the rows in each table.
    """
    return _sqlite.raw_tables(
        sqlite_callback = chembl_sqlite,
        sqlite = False,
        return_df = False
    )
