#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#
#  This file is part of the `bioigraph` python module
#
#  Copyright (c) 2014-2015 - EMBL-EBI
#
#  File author(s): Dénes Türei (denes@ebi.ac.uk)
#
#  Distributed under the GNU GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://www.ebi.ac.uk/~denes
#

import json
import sys
import bs4

import progress
import dataio

class Unichem(object):
    
    def __init__(self):
        sys.stdout.write('\n\tType `Unichem_instance.usage()` to get help.\n\n')
        sys.stdout.flush()
        # from unichem id to db name
        self.uc_dict = {
            '1': 'chembl',
            '2': 'drugbank',
            '3': 'pdb',
            '4': 'iuphar',
            '5': 'pubchem_dotf',
            '6': 'kegg_ligand',
            '7': 'chebi',
            '8': 'nih_ncc',
            '9': 'zinc',
            '10': 'emolecules',
            '11': 'ibm',
            '12': 'atlas',
            '13': 'patents',
            '14': 'fdasrs',
            '15': 'surechembl',
            '17': 'pharmgkb',
            '18': 'hmdb',
            '20': 'selleck',
            '21': 'pubchem_tpharma',
            '22': 'pubchem',
            '23': 'mcule',
            '24': 'nmrshiftdb2',
            '25': 'lincs'
        }
        # from db name to unichem id
        self.name_dict = {}
        for k,v in self.uc_dict.iteritems():
            self.name_dict[v] = k
        self.url_stem = 'https://www.ebi.ac.uk/unichem/rest/src_compound_id/'
        self.chembl_url = 'http://www.ebi.ac.uk/chemblws/compounds/smiles/{0}.json'
        self.cpd_search = 'http://www.ebi.ac.uk/unichem/rest/cpd_search/{0}/{1}/{2}'
        self.result = {}
    
    def usage(self):
        msg = '''
        List of identifier types can be read above. 
        To query UniChem, give either names or numbers of the
        ID types you wish to translate from and to.
        E.g. 
        >>> u = unichem.Unichem()
        >>> u.translate('pubchem','chembl',list_of_pubchems)
        
        Other option to search is connectivity search from UniChem.
        A-G parameters can be defined optionally. See description at
        https://www.ebi.ac.uk/unichem/info/widesearchInfo
        >>> u.connectivity_search(list_of_zincs,'zinc',parameters=[1,0,0,0,0,1,0])
        
        It's possible to get ChEMBL IDs from SMILEs via ChEMBL web service:
        >>> u.smiles2chembl(list_of_smiles)
        
        Find the dict in `u.result`. Untranslated items have value `None`.
        Every call overwrites previous result!
        '''
        sys.stdout.write('\n')
        for k,v in self.uc_dict.iteritems():
            sys.stdout.write('\t%s\t%s\n' % (k,v))
        sys.stdout.write(msg+'\n')
        sys.stdout.flush()
    
    def translate(self,source,target,lst):
        self.result = {}
        source = str(source) if type(source) is int else self.name_dict[source]
        target = str(target) if type(target) is int else self.name_dict[target]
        prg = progress.Progress(total=len(lst),name='Translating compound identifiers',
                                interval=1)
        for comp in lst:
            url = '/'.join([self.url_stem,comp,source,target])
            result = dataio.curl(url)
            self.result[comp] = []
            if result is not None:
                data = json.loads(result)
                for d in data:
                    self.result[comp].append(d['src_compound_id'])
            prg.step()
        prg.terminate()
    
    def smiles2chembl(self,smiles):
        self.result = {}
        prg = progress.Progress(total=len(smiles),name='Translating SMILEs',
                                interval=1)
        for sml in smiles:
            url = self.chembl_url.format(sml)
            result = dataio.curl(url)
            self.result[sml] = []
            if result is not None:
                try:
                    data = json.loads(result)
                    for d in data['compounds']:
                            this_smile = d['smiles']
                            this_chembl = d['chemblId']
                            #if this_smile == sml:
                            self.result[sml].append(this_chembl)
                except ValueError:
                    soup = bs4.BeautifulSoup(result)
                    compounds = soup.find_all('compound')
                    if compounds is not None:
                        for compound in compounds:
                            this_smile = compound.find('smiles').text
                            this_chembl = compound.find('chemblid').text
                            #if this_smile == sml:
                            self.result[sml].append(this_chembl)
            prg.step()
        prg.terminate()
    
    def connectivity_search(self,id_list,id_type,parameters=[1,0,0,0,0,1,0]):
        '''
        [1,0,0,0,0,1,0,  1]
        '''
        '''
        parameters is a list of parameters A-H as described in 
        https://www.ebi.ac.uk/unichem/info/widesearchInfo
        '''
        parameters.append(1) # H parameter must be 1 to process the result
        parameters = [str(i) for i in parameters]
        self.result = {}
        id_type = str(id_type) if type(id_type) is int else self.name_dict[id_type]
        prg = progress.Progress(total=len(id_list),name='Connectivity search',
                                interval=1)
        for i in id_list:
            prg.step()
            url = self.cpd_search.format(i,id_type,'/'.join(parameters))
            result = dataio.curl(url)
            self.result[i] = []
            if result is not None:
                data = json.loads(result)
                for k,v in data.iteritems():
                    for j in range(1,len(v)):
                        self.result[i].append(v[j][0])
            self.result[i] = list(set(self.result[i]))
        prg.terminate()