#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright (c) 2014-2016 - EMBL-EBI
#
#  File author(s): Dénes Türei (denes@ebi.ac.uk)
#
#  Distributed under the GNU GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://www.ebi.ac.uk/~denes
#

# mapping input methods

import pypath.uniprot_input as uniprot_input
import pypath.curl as curl
import pypath.urls as urls


def get_uniprot_sec(organism=9606):
    if organism is not None:
        proteome = uniprot_input.all_uniprots(organism=organism)
        proteome = set(proteome)
    sec_pri = []
    url = urls.urls['uniprot_sec']['url']
    c = curl.Curl(url, silent=False, large=True)
    data = c.result
    return filter(
        lambda line: len(line) == 2 and (organism is None or line[1] in proteome),
        map(lambda i: i[1].decode('utf-8').split(),
            filter(lambda i: i[0] >= 30, enumerate(data))))
