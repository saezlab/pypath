#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2019
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#                  Nicolàs Palacio
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

import os

import pypath.settings as settings


def get_cachedir(cachedir = None):
    """
    Ensures the cache directory exists and returs its path.
    """
    
    cachedir = cachedir or settings.get('cachedir')
    
    os.makedirs(cachedir, exist_ok = True)
    
    return cachedir
