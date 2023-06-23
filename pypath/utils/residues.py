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
try:
    import urllib2
except:
    import urllib.request as urllib2

import gzip
import bs4
try:
    from cStringIO import StringIO
except:
    try:
        from StringIO import StringIO
    except:
        from io import StringIO


class ResidueMapper(object):
    """
    This class stores and serves the PDB --> UniProt 
    residue level mapping. Attempts to download the 
    mapping, and stores it for further use. Converts 
    PDB residue numbers to the corresponding UniProt ones.
    """

    def __init__(self):
        self.url = 'http://pdb.org/pdb/rest/das/pdb_uniprot_mapping/alignment?query=%s'
        self.pdb_lst = 'ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/'\
            'pdb_chain_uniprot.tsv.gz'
        self.uniprot_pdb = None
        self.clean()
        self.download_errors = []

    def load_mapping(self, pdb):
        data = None
        non_digit = re.compile(r'[^\d.-]+')
        pdb = pdb.lower()
        url = self.url % pdb
        for i in range(5):
            try:
                data = urllib2.urlopen(url, timeout=60)
                break
            except:
                continue
        if not data:
            self.download_errors.append(pdb)
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
            if uniprot not in mapper:
                mapper[uniprot] = {}
            if chain not in mapper[uniprot]:
                mapper[uniprot][chain] = {}
            mapper[uniprot][chain][uniprotend] = {
                'pdbstart': pdbstart,
                'pdbend': pdbend,
                'uniprotstart': uniprotstart
            }
        self.mappers[pdb] = mapper

    def chains(self, chains):
        if type(chains) in [str, unicode]:
            chains = [chains]
        if type(chains) is list:
            chains = list(set(chains))
        return chains

    def pdb2uniprot(self, pdb, resnum, chains=None):
        chains = self.chains(chains)
        results = {}
        pdb = pdb.lower()
        if pdb not in self.mappers:
            self.load_mapping(pdb)
        if pdb in self.mappers:
            for ch, data in self.mappers[pdb].iteritems():
                if len(ch) == 1 and (chains is None or ch in chains):
                    pdbends = data.keys()
                    if resnum <= max(pdbends):
                        pdbend = min([
                            x for x in [e - resnum for e in pdbends] if x >= 0
                        ]) + resnum
                        seg = data[pdbend]
                        if seg['pdbstart'] <= resnum:
                            offset = seg['uniprotstart'] - seg['pdbstart']
                            residue = {
                                'resnum': resnum + offset,
                                'offset': offset,
                                'uniprot': seg['uniprot']
                            }
                            results[ch] = residue
        return results

    def uniprot2pdb(self, uniprot, resnum, chains=None, pdbs=None):
        chains = self.chains(chains)
        if self.uniprot_pdb is None:
            self.get_pdb_chains()
        results = {}
        # one uniprot can occure in more pdbs, first
        # we need to find out, which pdb files should we look at:
        if pdbs is None:
            pdbs = []
            if uniprot in self.uniprot_pdb:
                for updb in self.uniprot_pdb[uniprot]:
                    pdbs.append(updb['pdb'])
        elif type(pdbs) in [str, unicode]:
            pdbs = [pdbs]
        pdbs = list(set(pdbs))
        # now find the residue number in each of the pdb's:
        for pdb in pdbs:
            if pdb not in self.mappers:
                self.load_mapping(pdb)
            if pdb in self.mappers and uniprot in self.mappers[pdb]:
                for ch, up in self.mappers[pdb][uniprot].iteritems():
                    if chains is None or ch in chains:
                        uniprotends = up.keys()
                        if resnum <= max(uniprotends):
                            uniprotend = min([
                                x for x in [e - resnum for e in uniprotends]
                                if x >= 0
                            ]) + resnum
                            seg = up[uniprotend]
                            if seg['uniprotstart'] <= resnum:
                                offset = seg['pdbstart'] - seg['uniprotstart']
                                residue = {
                                    'resnum': resnum + offset,
                                    'offset': offset
                                }
                                if pdb not in results:
                                    results[pdb] = {}
                                results[pdb][ch] = residue
        return results

    def get_residue(self, ac, resnum, chains=None, pdbs=None):
        if len(ac.strip()) == 4:
            return self.pdb2uniprot(ac, resnum, chains)
        else:
            return self.uniprot2pdb(ac, resnum, chains, pdbs)

    def clean(self):
        '''
        Removes cached mappings, freeing up memory.
        '''
        self.mappers = {}
        self.uniprot_pdb = None
        self.pdb_uniprot = None

    def get_pdb_chains(self):
        gzfile = urllib2.urlopen(self.pdb_lst)
        buff = StringIO(gzfile.read())
        chains = gzip.GzipFile(fileobj=buff, mode='rb').read()
        chains = chains.replace('\r', '').split('\n')
        del chains[0]
        del chains[0]
        self.pdb_uniprot = {}
        self.uniprot_pdb = {}
        non_digit = re.compile(r'[^\d.-]+')
        for l in chains:
            l = l.split('\t')
            if len(l) > 8:
                if l[0] not in self.pdb_uniprot:
                    self.pdb_uniprot[l[0]] = {}
                self.pdb_uniprot[l[0]][l[1]] = {
                    'uniprot': l[2],
                    'chain_beg': int(non_digit.sub('', l[3])),
                    'chain_end': int(non_digit.sub('', l[4])),
                    'pdb_beg': int(non_digit.sub('', l[5])),
                    'pdb_end': int(non_digit.sub('', l[6])),
                    'uniprot_beg': int(non_digit.sub('', l[7])),
                    'uniprot_end': int(non_digit.sub('', l[8]))
                }
                if self.pdb_uniprot[l[0]][l[1]]['pdb_end'] - \
                        self.pdb_uniprot[l[0]][l[1]]['pdb_beg'] == \
                        self.pdb_uniprot[l[0]][l[1]]['uniprot_end'] - \
                        self.pdb_uniprot[l[0]][l[1]]['uniprot_beg']:
                    self.pdb_uniprot[l[0]][l[1]]['offset'] = \
                        (self.pdb_uniprot[l[0]][l[1]]['uniprot_beg'] -
                         self.pdb_uniprot[l[0]][l[1]]['pdb_beg'])
                else:
                    self.pdb_uniprot[l[0]][l[1]]['offset'] = None
                if l[2] not in self.uniprot_pdb:
                    self.uniprot_pdb[l[2]] = []
                self.uniprot_pdb[l[2]].append({
                    'pdb': l[0],
                    'chain': l[1],
                    'chain_beg': int(non_digit.sub('', l[3])),
                    'chain_end': int(non_digit.sub('', l[4])),
                    'pdb_beg': int(non_digit.sub('', l[5])),
                    'pdb_end': int(non_digit.sub('', l[6])),
                    'uniprot_beg': int(non_digit.sub('', l[7])),
                    'uniprot_end': int(non_digit.sub('', l[8])),
                    'offset': self.pdb_uniprot[l[0]][l[1]]['offset']
                })
