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

import download_manager as dm
from pypath.resources import urls


def download(url, fmt = 'open', **kwargs):
    """

    Args:
        fmt:
            "open" or "item" or "path".
    """

    res, dataset = url.split('/')
    url = urls.urls.get(res, {}).get(dataset, None) or url

    desc, item, downloader, path = _get_manager()._download(url, **kwargs)

    if fmt == 'open':

        return item.open(large = True)

    elif fmt == 'path':

        return path

    elif fmt == 'item':

        return item


def _get_manager():

    if 'DLMANAGER' not in globals():

        globals()['DLMANAGER'] = dm.DownloadManager(pkg = 'pypath')

    return globals()['DLMANAGER']
