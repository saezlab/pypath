#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright (c) 2014-2017 - EMBL-EBI
#
#  File author(s): Dénes Türei (denes@ebi.ac.uk)
#
#  Distributed under the GNU GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://www.ebi.ac.uk/~denes
#

from future.utils import iteritems

import re

import pypath.urls as urls
import pypath.curl as curl


def all_uniprots(organism=9606, swissprot=None):
    swissprot = 'YES' if swissprot == True else swissprot
    rev = '' if not swissprot else ' AND reviewed: %s' % swissprot
    url = urls.urls['uniprot_basic']['url']
    post = {
        'query': 'organism:%s%s' % (str(organism), rev),
        'format': 'tab',
        'columns': 'id'
    }
    c = curl.Curl(url, post=post, silent=False)
    data = c.result
    return list(
        filter(lambda x: len(x) > 0,
               map(lambda l: l.strip(), data.split('\n')[1:])))
