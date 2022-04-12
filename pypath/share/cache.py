#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#  (Planned for) centrally handling cache for all databases/resources.
#
#  Copyright
#  2014-2022
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

import os

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

    Args:
        key (str): A label for a named cache item. These are typically
            processed data dumped by the processing function for a quicker
            loading at next use.
    """

    if key in settings.settings.in_cachedir:

        return os.path.join(get_cachedir(), settings.get(key))
