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

import re

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.inputs.common as inputs_common
import pypath.formats.obo as obo


def get_pro():
    """
    Downloads and preprocesses the Protein Ontology.
    """
    
    url = urls.urls['pro']['url']
    
    reader = obo.Obo(url)
    
    return reader


def pro_mapping(target_id_type = 'UniProtKB', uniprot_isoforms = False):
    
    reid = re.compile(r'^(?:([\w_]+):)?([\w_-]+)')
    
    result = []
    
    url = urls.urls['pro']['mapping']
    
    c = curl.Curl(url, silent = False, large = True)
    
    _ = next(c.result)
    
    target_id_type = (
        'UniProtKB'
            if target_id_type == 'uniprot' else
        target_id_type
    )
    
    for line in c.result:
        
        line = line.split('\t')
        
        pro_id = line[0]
        
        id_type, target_id = reid.match(line[1]).groups()
        
        id_type = id_type or target_id.split('_')[0]
        
        if id_type == target_id_type:
            
            if uniprot_isoforms and target_id_type == 'UniProtKB':
                
                target_id, isoform = inputs_common._try_isoform(target_id)
            
            result.append((pro_id, target_id))
    
    return result
