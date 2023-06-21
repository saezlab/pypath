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
import pypath.share.curl as curl
import pypath.utils.mapping as mapping


def cellcellinteractions_annotations():
    
    
    CellcellinteractionsAnnotation = collections.namedtuple(
        'CellcellinteractionsAnnotation',
        [
            'mainclass',
        ]
    )
    
    
    url = urls.urls['cellcellinteractions']['url']
    
    c = curl.Curl(url, silent = False, large = True)
    
    _ = next(c.result)
    
    result = collections.defaultdict(set)
    
    for row in c.result:
        
        row = row.strip('\r\n').split('\t')
        
        uniprots = mapping.map_name(row[0], 'genesymbol', 'uniprot')
        classes = row[1].split('/')
        
        for uniprot in uniprots:
            
            for cls in classes:
                
                result[uniprot].add(
                    CellcellinteractionsAnnotation(mainclass = cls)
                )
    
    return dict(result)
