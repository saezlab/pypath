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
    Uses GitHub CSV file as alternative source since IPI database is discontinued.
    """
    
    result = collections.defaultdict(set)
    
    # Use GitHub CSV file as alternative source
    url = 'https://raw.githubusercontent.com/sacdallago/IPI_to_UniProt/master/ipi_uniprot_mapping.csv'
    
    c = curl.Curl(url, large = True, silent = False)
    
    if c.result is None:
        return dict(result)
    
    # Skip header row
    rows = iter(c.result)
    next(rows, None)
    
    for row in rows:
        
        row = row.strip('\n\r').split(',')
        
        if len(row) < 2:
            continue
        
        ipi_id = row[0].strip()
        uniprot = row[1].strip()
        
        if ipi_id and uniprot:
            result[ipi_id].add(uniprot)
    
    return dict(result)


def _ipi_uniprot_pairs(*args, **kwargs):
    
    for ipi, uniprots in iteritems(ipi_uniprot()):
        
        for uniprot in uniprots:
            
            yield ipi, uniprot
