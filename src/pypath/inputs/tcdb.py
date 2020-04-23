#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2020
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#                  Nicolàs Palacio
#                  Olga Ivanova
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

from future.utils import iteritems

import re
import collections
import bs4

import pypath.inputs.common as inputs_common
import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.utils.mapping as mapping


def tcdb_families():
    
    retag = re.compile(r'<.*>')
    
    url = urls.urls['tcdb']['url_families']
    
    c = curl.Curl(url, large = False, silent = False)
    
    lines = bs4.BeautifulSoup(c.result, features = 'lxml').find('p').text
    
    return dict(
        (
            tcid,
            family.replace('\t', ' ')
        )
        for tcid, family in
        (
            retag.sub('', line.strip()).split('\t', maxsplit = 1)
            for line in lines.strip().split('\n')
        )
    )
