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

import pypath.resources.urls as urls
import pypath.share.curl as curl


def depod_interactions(organism = 9606):

    url = urls.urls['depod']['urls'][1]
    c = curl.Curl(url, silent = False, large = True, encoding = 'iso-8859-1')
    data = c.result
    result = []
    i = []
    lnum = 0

    for l in data:

        if lnum == 0:
            lnum += 1
            continue
        l = l.replace('\n', '').replace('\r', '')
        l = l.split('\t')
        specA = int(l[9].split(':')[1].split('(')[0])
        specB = int(l[10].split(':')[1].split('(')[0])
        
        if organism is None or (specA == organism and specB == organism):
            
            pm = l[8].replace('pubmed:', '')
            sc = l[14].replace('curator score:', '')
            ty = l[11].split('(')[1].replace(')', '')
            l = [l[0], l[1]]
            interaction = ()
            
            for ll in l:
                
                ll = ll.split('|')
                uniprot = ''
                for lll in ll:
                    nm = lll.split(':')
                    u = nm[1].strip()
                    if nm[0] == 'uniprotkb' and len(u) == 6:
                        uniprot = u
                interaction += (uniprot, )
            
            interaction += (pm, sc, ty)
            if len(interaction[0]) > 1 and len(interaction[1]) > 1:
                i.append(interaction)
        
        lnum += 1

    return i
