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
import bs4

try:
    import urllib2
except:
    import urllib.request as urllib2

import pypath.resources.data_formats as data_formats


class ResidueMapper(object):
    """
    This class stores and serves the PDB --> UniProt 
    residue level mapping. Attempts to download the 
    mapping, and stores it for further use. Converts 
    PDB residue numbers to the corresponding UniProt ones.
    """


    def __init__(self):
        
        self.clean()


    def load_mapping(self, pdb):
        
        non_digit = re.compile(r'[^\d.-]+')
        pdb = pdb.lower()
        url = data_formats.urls['pdb_align']['url'] + pdb
        data = urllib2.urlopen(url)
        mapper = {}
        soup = bs4.BeautifulSoup(data.read())
        
        for block in soup.find_all('block'):
            
            seg = block.find_all('segment')
            chain = seg[0]['intobjectid'].split('.')[1]
            uniprot = seg[1]['intobjectid']
            pdbstart = int(non_digit.sub('', seg[0]['start']))
            pdbend = int(non_digit.sub('', seg[0]['end']))
            uniprotstart = int(non_digit.sub('', seg[1]['start']))
            uniprotend = int(non_digit.sub('', seg[1]['end']))
            if chain not in mapper:
                mapper[chain] = {}
            mapper[chain][pdbend] = {
                'uniprot': uniprot,
                'pdbstart': pdbstart,
                'uniprotstart': uniprotstart,
                'uniprotend': uniprotend
            }
        
        self.mappers[pdb] = mapper


    def get_residue(self, pdb, resnum, chain = None):
        
        pdb = pdb.lower()
        
        if pdb not in self.mappers:
            self.load_mapping(pdb)
            
        if pdb in self.mappers:
            
            for chain, data in self.mappers[pdb].iteritems():
                
                pdbends = data.keys()
                
                if resnum <= max(pdbends):
                    
                    pdbend = min(
                        [x for x in [e - resnum for e in pdbends]
                         if x >= 0]) + resnum
                    seg = data[pdbend]
                    
                    if seg['pdbstart'] <= resnum:
                        
                        offset = seg['uniprotstart'] - seg['pdbstart']
                        residue = {
                            'resnum': resnum + offset,
                            'offset': offset,
                            'uniprot': seg['uniprot'],
                            'chain': chain
                        }
                        
                        return residue
        
        return None


    def clean(self):
        """
        Removes cached mappings, freeing up memory.
        """
        
        self.mappers = {}
