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
        self.url_stem = 'https://www.ebi.ac.uk/unichem/rest/src_compound_id'
        self.inchi_stem = 'https://www.ebi.ac.uk/unichem/rest/inchikey/%s'
        self.chembl_url = 'http://www.ebi.ac.uk/chemblws/compounds/smiles/{0}.json'
        self.cpd_search = 'http://www.ebi.ac.uk/unichem/rest/{0}/{1}/{2}{3}'
        self.result = {}
    
    def usage(self):
        msg = '''
        List of identifier types can be read above. 
        To query UniChem, give either names or numbers of the
        ID types you wish to translate from and to.
        E.g. 
        >>> u = unichem.Unichem()
        >>> u.translate('pubchem', 'chembl', list_of_pubchems)
        
        Additional ways of translation are from SMILEs to ChEMBL IDs, and
        from InChiKeys to any ID. These are translated not by the UniChem,
        but by the ChEMBL webservice.
        >>> u.translate('smiles', 'chembl', list_of_smiles)
        >>> u.translate('inchikey', 'chembl', list_of_inchikeys)
        
        Other option to search is connectivity search from UniChem.
        A-G parameters can be defined optionally. See description at
        https://www.ebi.ac.uk/unichem/info/widesearchInfo
        >>> u.connectivity_search(list_of_zincs, 'zinc', parameters=[1,0,0,0,0,1,0])
        
        InChiKeys can be used in connectivity search too:
        >>> u.connectivity_search(list_of_inchikeys, 'inchikey', parameters=[1,0,0,0,0,1,0])
        
        You can also call directly functions accessing ChEMBL webservice, with the 
        same result as you would call `translate()` or `connectivity_search()`:
        >>> u.smiles2chembl(list_of_smiles)
        >>> u.inchikey2anything('chembl', list_of_inchikeys)
        
        Find the dict in `u.result`. Untranslated items have value `None`.
        Every call overwrites previous result!
        '''
        sys.stdout.write('\n')
        for k,v in self.uc_dict.iteritems():
            sys.stdout.write('\t%s\t%s\n' % (k,v))
        sys.stdout.write(msg+'\n')
        sys.stdout.flush()
    
    def translate(self,source,target,lst):
        if source == 'inchikey':
            self.inchikey2anything(target, lst)
            return None
        if source == 'smiles':
            self.smiles2chembl(lst)
            return None
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
    
    def inchikey2anything(self, target, lst):
        self.result = {}
        target = str(target) if type(target) is int else self.name_dict[target]
        prg = progress.Progress(total=len(lst),name='Translating InChi-Keys',
            interval=1)
        for inchik in lst:
            url = self.inchi_stem%inchik
            result = dataio.curl(url)
            if result is not None:
                data = json.loads(result)
                self.result[inchik] = [d['src_compound_id'] \
                    for d in data if d['src_id'] == target]
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
    
    def connectivity_search(self, id_list, id_type, parameters=[1,0,0,0,0,1,0]):
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
        if id_type == 'inchikey':
            id_type = ''
            method = 'key_search'
        elif id_type == 'smiles':
            self.result = None
            return None
        else:
            id_type = str(id_type) if type(id_type) is int else self.name_dict[id_type]
            id_type = '%s/'%id_type
            method = 'cpd_search'
        prg = progress.Progress(total=len(id_list),name='Connectivity search',
                                interval=1)
        for i in id_list:
            prg.step()
            url = self.cpd_search.format(method, i, id_type, '/'.join(parameters))
            result = dataio.curl(url)
            self.result[i] = []
            if result is not None:
                data = json.loads(result)
                for k,v in data.iteritems():
                    for j in range(1,len(v)):
                        self.result[i].append(v[j][0])
            self.result[i] = list(set(self.result[i]))
        prg.terminate()