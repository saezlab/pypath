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

from __future__ import annotations

from typing import Any

import os
import hashlib
import pathlib as pl

import pypath_common._misc as _common
import pypath.share.settings as settings


def get_cachedir(cachedir: str | pl.Path | None = None) -> pl.Path:
    """
    Ensures the cache directory exists and returns its path.
    """

    cachedir = pl.Path(
        settings.get(
            'cachedir',
            override = cachedir,
            default = settings.settings._user_cache_dir,
        ),
    )

    os.makedirs(cachedir, exist_ok = True)

    return cachedir


def cache_item(key):
    """
    For a key of a cache item returns its path. It does not mean the file
    actually exists.

    Args
        key (str): A label for a named cache item. These are typically
            processed data dumped by the processing function for a quicker
            loading at next use.
    """

    if key in settings.settings.in_cachedir:

        return os.path.join(get_cachedir(), settings.get(key))


def cache_path(key: Any) -> str:
    """
    Path to a cache item identified by a unique key.

    Similar to `cache_item`, but instead of using a key registered in the
    module config, it processes the key into an MD5 hash. The key can be
    anything: a string, a tuple, etc.
    """

    return os.path.join(
        get_cachedir(),
        hashlib.md5(str(key).encode('utf-8')).hexdigest(),
    )
