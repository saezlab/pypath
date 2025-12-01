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
import shutil
import os

import pypath.share.cache as cache
from pypath.share import curl, cache
from pypath.formats import sqlite as _sqlite
from pypath.resources import urls

from ._common import _log, _show_tables

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

    path_in_tar = (
        f'chembl_{version:02d}/chembl_{version:02d}_sqlite/chembl_{version:02d}.db'
        )
    url = urls.urls['chembl']['sqlite'] % (version, version)

    def _download() -> curl.Curl:
        """Callback to download the tarballed ChEMBL database."""
        return curl.Curl(
            url,
            large=True,
            silent=False,
            files_needed=[path_in_tar],
            slow=True
        )

    def _extractor(curl_obj: curl.Curl) -> str:
        """Extractor to get the SQLite file from the tarball."""
        return curl_obj.result[path_in_tar]

    return _sqlite.download_sqlite(
        download_callback=_download,
        extractor=_extractor,
        database = 'ChEMBL',
        version = f'{version:02}',
        connect = connect,
    )
