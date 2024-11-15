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

import pypath.share.curl as curl
import pypath.resources.urls as urls


def abs_interactions():
    """
    TF-target (transcriptional regulation) interactions from the ABS
    database.
    """

    result = []
    url = urls.urls['abs']['url']

    dmanager = dm.DownloadManager(pkg='pypath')
    desc, item, downloader, dest = dmanager._download(url, backend='curl', ssl_verifypeer=0, ssl_verifyhost=0)
    print(item, dest)
    c = item.open()

    data = c.result
    data = [[x.replace('*', '') for x in xx.split('\t')]
            for xx in data.split('\n')]

    for d in data:

        if len(d) > 2:

            result.append([d[2], d[0]])

    return result
