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

import re
import collections

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.utils.taxonomy as taxonomy
import pypath.inputs.common as inputs_common


resite = re.compile('([A-Z])(\d+)')


IptmnetRecord = collections.namedtuple(
    'IptmnetRecord',
    [
        'enzyme',
        'substrate',
        'substrate_isoform',
        'ptm_type',
        'resaa',
        'resnum',
        'score',
    ]
)


def iptmnet(organism = 9606):
    
    ptm_url = urls.urls['iptmnet']['ptms']
    score_url = urls.urls['iptmnet']['scores']
    protein_url = urls.urls['iptmnet']['proteins']
    
    c = curl.Curl(protein_url, large = True, silent = False)
    
    pr_ids = {}
    
    for line in c.result:
        
        line = line.strip('\n\r').split('\t')
        
        if line[5]:
            
            pr_ids[line[5]] = line[0]
    
    c = curl.Curl(score_url, large = True, silent = False)
    
    scores = {}
    
    for line in c.result:
        
        line = line.strip('\n\r').split('\t')
        
        if not line[2]:
            
            continue
        
        site = resite.match(line[1])
        
        if not site:
            
            continue
        
        resaa, resnum = site.groups()
        
        resnum = int(resnum)
        score = int(line[4])
        substrate, isoform = inputs_common._try_isoform(line[0])
        substrate = pr_ids[substrate] if substrate in pr_ids else substrate
        enzyme = pr_ids[line[2]] if line[2] in pr_ids else line[2]
        
        key = (
            enzyme,
            substrate,
            isoform,
            line[3].lower(), # PTM type
            resaa,
            resnum,
        )
        
        scores[key] = score
    
    c = curl.Curl(ptm_url, large = True, silent = False)
    
    for line in c.result:
        
        line = line.strip('\n\r').split('\t')
        
        if not line or not line[6]:
            
            continue
        
        ncbi_tax_id = taxonomy.ensure_ncbi_tax_id(line[4])
        
        if organism and ncbi_tax_id != organism:
            
            continue
        
        substrate, isoform = inputs_common._try_isoform(line[2])
        substrate = pr_ids[substrate] if substrate in pr_ids else substrate
        ptm_type = line[0].lower()
        enzyme = pr_ids[line[6]] if line[6] in pr_ids else line[6]
        refs = line[8].split(',')
        resnum, resaa = resite.match(line[5]).groups()
        
        key = (
            enzyme,
            substrate,
            isoform,
            ptm_type,
            resaa,
            resnum,
        )
        
        score = scores[key] if key in scores else None
        
        yield IptmnetRecord(
            enzyme = enzyme,
            substrate = substrate,
            substrate_isoform = isoform,
            ptm_type = ptm_type,
            resaa = resaa,
            resnum = resnum,
            score = score,
        )
