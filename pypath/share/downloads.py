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

"""
Shared DownloadManager instance for pypath.
"""

from pathlib import Path
from download_manager import DownloadManager


# Pypath data directory in the workspace
DATA_DIR = Path(__file__).parent.parent.parent.parent / 'pypath-data'


def get_download_manager() -> DownloadManager:
    """
    Get the shared DownloadManager instance configured with pypath's data folder.

    Returns:
        DownloadManager: Configured download manager instance.
    """
    return DownloadManager(data_folder=str(DATA_DIR))


# Singleton instance
dm = get_download_manager()
