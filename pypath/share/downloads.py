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

import os
from pathlib import Path
from typing import Optional, List

from dotenv import load_dotenv

import pypath.share.settings as settings
from dlmachine import DownloadManager
from cachedir._open import Opener


def _load_env_file() -> None:
    """Load environment variables from the nearest local `.env` file."""

    for base in (Path.cwd(), *Path.cwd().parents):

        env_path = base / '.env'

        if env_path.exists():
            load_dotenv(env_path, override=False)
            break


def _resolve_data_dir() -> Path:
    """Resolve download storage directory for inputs_v2 raw files."""

    _load_env_file()

    configured = (
        os.environ.get('PYPATH_DOWNLOAD_DATADIR') or
        settings.get('cachedir')
    )
    data_dir = Path(configured)
    data_dir.mkdir(parents=True, exist_ok=True)

    return data_dir


# Kept as module-level constant for compatibility with modules importing DATA_DIR
DATA_DIR = _resolve_data_dir()


def get_download_manager() -> DownloadManager:
    """
    Get the shared DownloadManager instance configured with pypath's data folder.

    Returns:
        DownloadManager: Configured download manager instance.
    """
    return DownloadManager(path=str(_resolve_data_dir()), config={'backend': 'requests'})


# Singleton instance
dm = get_download_manager()


def download_and_open(
        url: str,
        filename: str,
        subfolder: str,
        large: bool = True,
        encoding: str = 'utf-8',
        default_mode: str = 'r',
        ext: Optional[str] = None,
        needed: Optional[List[str]] = None,
        **download_kwargs,
    ) -> Opener:
    """
    Download a file and open it with curl.Curl-compatible interface.

    This function combines download with automatic extraction/opening,
    returning an Opener instance that mimics curl.Curl behavior.

    Args:
        url: URL to download
        filename: Filename to save as
        subfolder: Subfolder in the data directory
        large: If True, return file handles for streaming (default: True)
        encoding: Text encoding (default: 'utf-8')
        default_mode: File mode 'r' for text, 'rb' for binary (default: 'r')
        ext: File extension for compression detection ('zip', 'gz', 'tar.gz', etc.)
        needed: For archives, list of specific files to extract (default: all)
        **download_kwargs: Additional arguments passed to dm.download() (e.g., query, post)

    Returns:
        Opener instance with a .result attribute containing:
        - For zip/tar files: dict mapping filenames to file handles (if large=True)
        - For gz files: file handle
        - For plain files: file handle

    Example:
        >>> # For zip files - returns dict like curl.Curl
        >>> opener = download_and_open(url, 'data.zip', 'mydata')
        >>> for filename, handle in opener.result.items():
        >>>     # process file

        >>> # For plain files
        >>> opener = download_and_open(url, 'data.txt', 'mydata')
        >>> for line in opener.result:
        >>>     # process line
    """

    data_dir = _resolve_data_dir()

    # Download the file to a deterministic path
    target_dir = data_dir / subfolder
    target_dir.mkdir(parents=True, exist_ok=True)
    file_path = target_dir / filename
    dm.download(url, dest=str(file_path), **download_kwargs)

    # Use Opener to handle extraction/opening
    # Return the opener itself so it stays alive and keeps files open
    opener = Opener(
        path=str(file_path),
        ext=ext,
        needed=needed,
        large=large,
        default_mode=default_mode,
        encoding=encoding,
    )

    return opener
