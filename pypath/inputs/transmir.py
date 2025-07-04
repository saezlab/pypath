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

import collections

import pypath.resources.urls as urls
from download_manager import DownloadManager



def transmir_interactions():
    """
    TF-miRNA gene intractions from TransmiR.
    """

    url = urls.urls['transmir']['url']
    
    # Initialize download manager with cache
    dm = DownloadManager(pkg='pypath')
    
    # Use _download to get the cache item
    desc, item, downloader, path = dm._download(
        url,
        dest=True,  # Download to cache/file
    )

    TransmirInteraction = collections.namedtuple(
        'TransmirInteraction',
        [
            'tf_genesymbol',
            'mirna',
            'effect',
            'pubmed',
        ]
    )

    result = []
    
    # Use the cache item's open method which handles decompression automatically
    if item and item.status == 3:  # Status.READY = 3
        opener = item.open(large=True, encoding='iso-8859-1')
        
        if opener and opener.result:
            # opener.result is an iterator when large=True
            for l in opener.result:
                l = l.strip()
                if l:  # Skip empty lines
                    l = l.split('\t')
                    
                    if len(l) >= 6:  # Ensure we have enough fields
                        result.append(
                            TransmirInteraction(
                                tf_genesymbol = l[0],
                                mirna = l[1],
                                effect = l[4].split('(')[0],
                                pubmed = l[5],
                            )
                        )

    return result
