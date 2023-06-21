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
import pypath.utils.taxonomy as taxonomy


def gpcrdb_annotations(organism = 9606):
    """
    :param int,str organism:
        Only human and mouse (9606 and 10090) are supported.
    """
    
    
    GpcrdbAnnotation = collections.namedtuple(
        'GpcrdbAnnotation',
        [
            'gpcr_class',
            'family',
            'subfamily',
        ]
    )
    
    organism = taxonomy.ensure_ncbi_tax_id(organism)
    
    if organism not in (9606, 10090):
        
        return {}
    
    i_uniprot = 31 if organism == 10090 else 15
    
    url = urls.urls['gpcrdb']['families']
    
    c = curl.Curl(url, silent = False, large = True)
    
    result = collections.defaultdict(set)
    
    for line in c.result:
        
        if line[0] != ' ':
            
            cls = line.split('|')[0].strip()
            family = None
            subfamily = None
        
        elif line[4] != ' ':
            
            family = line.strip()
            subfamily = None
            
        elif line[8] != ' ':
            
            subfamily = line.strip()
            
        else:
            
            line = line.strip().strip('"')
            
            if line.startswith('gpcr'):
                
                line = line.split('","')
                uniprot = line[i_uniprot]
                
                if uniprot:
                    
                    result[uniprot].add(
                        GpcrdbAnnotation(
                            gpcr_class = cls,
                            family = family,
                            subfamily = subfamily,
                        )
                    )
    
    return dict(result)
