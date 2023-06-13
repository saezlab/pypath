#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#  (Planned for) centrally handling cache for all databases/resources.
#
#  Copyright
#  2014-2023
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: Dénes Türei (turei.denes@gmail.com)
#           Nicolàs Palacio
#           Sebastian Lobentanzer
#           Erva Ulusoy
#           Olga Ivanova
#           Ahmet Rifaioglu
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

from typing import Any

import os
import hashlib

import pypath.share.settings as settings


def get_cachedir(cachedir = None):
    """
    Ensures the cache directory exists and returns its path.
    """

    cachedir = cachedir or settings.get('cachedir')

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
