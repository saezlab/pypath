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

from future.utils import iteritems

import collections

import pypath.inputs.common as inputs_common
import pypath.resources.urls as urls
import pypath.share.curl as curl


def ipi_uniprot():
    """
    Retrieves an IPI-UniProt mapping dictionary.
    """
    
    result = collections.defaultdict(set)
    
    url = urls.urls['ipi']['url']
    
    c = curl.Curl(url, large = True, silent = False)
    
    for row in c.result:
        
        row = row.strip('\n\r').split('\t')
        
        if len(row) < 3:
            
            continue
        
        ipi_id = row[2]
        
        uniprot, isoform = inputs_common._try_isoform(row[1])
        
        is_uniprot = (
            not any(
                uniprot.startswith(pref)
                for pref in ('NP_', 'OTTH', 'HIT', 'ENSP', 'XP_')
            )
        )
        
        if is_uniprot:
            
            result[ipi_id].add(uniprot)
    
    return dict(result)


def _ipi_uniprot_pairs(*args, **kwargs):
    
    for ipi, uniprots in iteritems(ipi_uniprot()):
        
        for uniprot in uniprots:
            
            yield ipi, uniprot
