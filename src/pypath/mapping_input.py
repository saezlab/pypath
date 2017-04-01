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

# mapping input methods

from future.utils import iteritems

import re

import pypath.uniprot_input as uniprot_input
import pypath.curl as curl
import pypath.urls as urls
import pypath.common as common


def get_uniprot_sec(organism=9606):
    """
    Downloads and processes the mapping between secondary and
    primary UniProt IDs.
    
    :param int organism: NCBI Taxonomy ID of the organism.
    """
    
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

def get_mirbase_aliases(organism = 9606):
    """
    Downloads and processes mapping tables from miRBase.
    """
    
    if type(organism) in common.charTypes:
        mborganism = organism
    elif organism not in common.mirbase_taxids:
        raise ValueError('Organism not known: %u. Try to pass miRBase '
                         'taxon prefix as string, e.g. `hsa`.' % organism)
    else:
        mborganism = common.mirbase_taxids[organism]
    
    mat = {}
    mir = {}
    
    url = urls.urls['mirbase']['aliases']
    c = curl.Curl(url, silent = False, large = True)
    
    for l in c.result:
        
        l = l.decode('utf-8').strip().strip(';').split('\t')
        
        if l[1][:3] != mborganism:
            continue
        
        d = mat if l[0][:5] == 'MIMAT' else mir
        
        if l[0] not in d:
            d[l[0]] = set([])
        
        for m in l[1].split(';'):
            d[l[0]].add(m)
    
    return mat, mir

def mirbase_mature(organism = 9606):
    
    mat, mir = get_mirbase_aliases(organism)
    
    result = {}
    
    for mimat, mmats in iteritems(mat):
        
        for mmat in mmats:
            
            yield mimat, mmat

def mirbase_precursor(organism = 9606):
    
    mat, mir = get_mirbase_aliases(organism)
    
    result = {}
    
    for mi, mpres in iteritems(mir):
        
        for mpre in mpres:
            
            yield mi, mpre

def mirbase_precursor_to_mature(organism = 9606):
    
    pre = mirbase_precursor(organism)
    ids = mirbase_ids(organism)
    
    _ids = {}
    
    for mpre, mmat in ids:
        
        if mpre not in _ids:
            _ids[mpre] = set([])
        
        _ids[mpre].add(mmat)
    
    result = {}
    
    for prename, mpres in pre:
        
        for mpre in mpres:
            
            if mpre in _ids:
                
                for mmat in _ids[mpre]:
                    
                    yield prename, mmat

def mirbase_ids(organism = 9606):
    
    reprename = re.compile(r'([-A-z]*[-]?\d+[a-z]*)(-\d*)')
    
    def get_pre_name(mat_name):
        return mat_name.replace(
            '*', '').replace(
            '-3p', '').replace(
            '-5p', '')
    
    mat, mir = get_mirbase_aliases(organism)
    
    mir = dict((k, set.union(set(reprename.sub(r'\1', vv) for vv in v), v))
               for k, v in iteritems(mir))
    
    mir = common.swap_dict(mir)
    
    mat = dict((k, set(get_pre_name(vv) for vv in v))
               for k, v in iteritems(mat))
    
    if (sum(sum(vv in mir for vv in v) for v in mat.values()) <
        sum(sum(vv.lower() in mir for vv in v) for v in mat.values())):
        
        mat = dict((k, set(vv.lower() for vv in v))
               for k, v in iteritems(mat))
    
    mat_mir = common.join_dicts(mat, mir)
    
    for ma, mis in iteritems(mat_mir):
        
        for mi in (mis if type(mis) not in common.simpleTypes else [mis]):
            
            yield ma, mi
