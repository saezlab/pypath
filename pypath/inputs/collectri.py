#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
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
#           Sophia Müller-Dott
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

import re
import collections

import pypath.resources.urls as urls
import pypath.share.curl as curl


def collectri_interactions(int_type = 'gene'):

    CollectriInteraction = collections.namedtuple(
        'CollectriInteraction',
        (
            'tf',
            'target',
            'effect',
            'TF_category',
            'resources',
            'pubmed',
        ),
    )

    url = urls.urls['collectri']['url']
    c = curl.Curl(
        url,
        silent = False,
        large = True,
    )

    result = []
    
    if int_type == 'miRNA':
        for l in c.result:
            l = l.strip().split(',')
            match = re.match('^MIR(\d+)$', l[1])
            
            if match:
                number = match.group(1)
                l[1] = 'hsa-miR-' + number
                result.append(
                    CollectriInteraction(l[0], l[1], int(l[2]), l[3], l[4], l[5])
                    )
                
        return result

    if int_type == 'gene':
        for l in c.result:
            l = l.strip().split(',')
            
            
            if l[0] != 'source':
                result.append(
                    CollectriInteraction(l[0], l[1], int(l[2]), l[3], l[4], l[5])
                    )
                
        return result