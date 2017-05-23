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

from future.utils import iteritems
from past.builtins import xrange, range, reduce

import os
import sys
import codecs
import re
import imp
import copy
import itertools

import urllib

if not hasattr(urllib, 'urlencode'):
    _urllib = urllib
    urllib = _urllib.parse

import json
try:
    import cPickle as pickle
except:
    import pickle

# from pypath:
import pypath.progress as progress
import pypath.logn as logn
import pypath.common as common
import pypath.maps as maps
import pypath.mysql as mysql
import pypath.urls as urls
import pypath.curl as curl
import pypath.mapping_input as mapping_input
import pypath.uniprot_input as uniprot_input
import pypath.input_formats as input_formats

__all__ = ['MappingTable', 'Mapper']

###
# functions to read and use mapping tables from UniProt, file, mysql or pickle
###


class MappingTable(object):

    def __init__(self,
                 one,
                 two,
                 typ,
                 source,
                 param,
                 ncbi_tax_id,
                 mysql=None,
                 log=None,
                 cache=False,
                 cachedir='cache',
                 uniprots = None):
        
        '''
        When initializing ID conversion tables for the first time
        data is downloaded from UniProt and read into dictionaries.
        It takes a couple of seconds. Data is saved to pickle 
        dumps, this way later the tables load much faster.
        '''
        
        self.param = param
        self.one = one
        self.two = two
        self.typ = typ
        self.maxlOne = None
        self.maxlTwo = None
        self.mysql = mysql
        self.cache = cache
        self.cachedir = cachedir
        self.mapping = {"to": {}, "from": {}}
        
        if log.__class__.__name__ != 'logw':
            self.session = common.gen_session_id()
            self.ownlog = logn.logw(self.session, 'INFO')
        else:
            self.ownlog = log
        
        if param is not None:
            
            self.mid = common.md5((one, two, self.param.bi, ncbi_tax_id))
            md5param = common.md5(json.dumps(self.param.__dict__))
            self.cachefile = os.path.join(self.cachedir, md5param)
            
            if self.cache and os.path.isfile(self.cachefile):
                self.mapping = pickle.load(open(self.cachefile, 'rb'))
            
            elif len(self.mapping['to']) == 0 or (
                    param.bi and len(self.mapping['from']) == 0):
                
                if os.path.exists(self.cachefile):
                    os.remove(self.cachefile)
                if source == "mysql":
                    self.read_mapping_mysql(param, ncbi_tax_id)
                elif source == "file":
                    self.read_mapping_file(param, ncbi_tax_id)
                elif source == "pickle":
                    self.read_mapping_pickle(param, ncbi_tax_id)
                elif source == "uniprot":
                    self.read_mapping_uniprot(param, ncbi_tax_id)
                elif source == "uniprotlist":
                    self.read_mapping_uniprot_list(param,
                                                   uniprots = uniprots,
                                                   ncbi_tax_id = ncbi_tax_id)
                
                if len(self.mapping['to']) and (
                        not param.bi or len(self.mapping['from'])):
                    pickle.dump(self.mapping, open(self.cachefile, 'wb'))

    def reload(self):
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)

    def cleanDict(self, mapping):
        for key, value in iteritems(mapping):
            mapping[key] = common.uniqList(value)
        return mapping

    def read_mapping_file(self, param, ncbi_tax_id = None):
        
        ncbi_tax_id = self.get_tax_id(ncbi_tax_id)
        
        if param.__class__.__name__ != "FileMapping":
            self.ownlog.msg(2, "Invalid parameter for read_mapping_file()",
                            'ERROR')
            return {}
        
        if (not os.path.exists(param.input) and
            not hasattr(mapping_input, param.input)):
            
            return {}
        
        if hasattr(mapping_input, param.input):
            
            toCall = getattr(mapping_input, param.input)
            inputArgs = param.inputArgs if hasattr(param, 'inputArgs') else {}
            infile = list(toCall(**inputArgs))
            
            total = sum([sys.getsizeof(i) for i in infile])
            
        else:
            infile = codecs.open(param.input, encoding='utf-8', mode='r')
            total = os.path.getsize(param.input)
        
        prg = progress.Progress(
            total=total, name="Reading from file", interval=18)
        
        lnum = 0
        lsum = 0
        mapping_o = {}
        mapping_i = {}
        
        for line in infile:
            
            if len(line) == 0:
                continue
            
            if lnum == 0 and param.header != 0:
                lnum += 1
                continue
            
            if type(line is list):
                prg.step(sys.getsizeof(line))
                
            else:
                line = line.decode('utf-8')
                prg.step(len(line))
                line = line.rstrip().split(param.separator)
            
            if len(line) > max([param.oneCol, param.twoCol]):
                
                if line[param.oneCol] not in mapping_o:
                    
                    mapping_o[line[param.oneCol]] = []
                mapping_o[line[param.oneCol]].append(line[param.twoCol])
                
                if param.bi:
                    
                    if line[param.twoCol] not in mapping_i:
                        
                        mapping_i[line[param.twoCol]] = []
                    mapping_i[line[param.twoCol]].append(line[param.oneCol])
            
            lnum += 1
        
        if hasattr(infile, 'close'):
            infile.close()
        
        self.mapping["to"] = mapping_o
        self.cleanDict(self.mapping["to"])
        
        if param.bi:
            self.mapping["from"] = mapping_i
            self.cleanDict(self.mapping["from"])
        
        prg.terminate()
    
    def read_mapping_uniprot_list(self, param, uniprots = None,
                                  ncbi_tax_id = None):
        
        mapping_o = {}
        mapping_i = {}
        
        ncbi_tax_id = param.ncbi_tax_id \
            if ncbi_tax_id is None else ncbi_tax_id
        
        if uniprots is None:
            uniprots = uniprot_input.all_uniprots(ncbi_tax_id,
                                                  swissprot = param.swissprot)
        
        if param.targetNameType != 'uniprot':
            utarget = self._read_mapping_uniprot_list('ACC',
                                                      param.target_ac_name,
                                                      uniprots)
            
            _ = utarget.readline()
            ac_list = list(map(lambda l:
                                   l.decode('ascii').split('\t')[1].strip(),
                                   utarget))
        else:
            ac_list = uniprots
        
        udata = self._read_mapping_uniprot_list(param.target_ac_name,
                                   param.ac_name,
                                   ac_list)
        
        _ = udata.readline()
        
        for l in udata:
            
            l = l.decode('ascii').strip().split('\t')
            
            if l[1] not in mapping_o:
                mapping_o[l[1]] = []
            
            mapping_o[l[1]].append(l[0])
            
            if param.bi:
                
                if l[0] not in mapping_i:
                    mapping_i[l[0]] = []
                
                mapping_i[l[0]].append(l[1])
        
        self.mapping["to"] = mapping_o
        self.cleanDict(self.mapping["to"])
        if param.bi:
            self.mapping["from"] = mapping_i
            self.cleanDict(self.mapping["from"])
    
    def _read_mapping_uniprot_list(self, source, target, ac_list):
        """
        Reads a mapping table from UniProt "upload lists" service.
        """
        
        url = urls.urls['uniprot_basic']['lists']
        post = {
            'from': source,
            'format': 'tab',
            'to': target,
            'uploadQuery': ' '.join(ac_list)
        }
        
        c = curl.Curl(url, post=post, large=True, silent = False)
        
        if c.result is None:
            for i in xrange(3):
                c = curl.Curl(url, post=post, large=True,
                              silent = False, cache = False)
                if c.result is not None:
                    break
            
            if c.result is None:
                sys.stdout.write('\t:: Error at downloading from UniProt.\n')
        
        return c.result
    
    def read_mapping_uniprot(self, param, ncbi_tax_id = None):
        """
        Downloads ID mappings directly from UniProt.
        See the names of possible identifiers here:
        http://www.uniprot.org/help/programmatic_access

        :param UniprotMapping param: UniprotMapping instance
        :param int ncbi_tax_id: Organism NCBI Taxonomy ID.
        """
        
        ncbi_tax_id = self.get_tax_id(ncbi_tax_id)
        resep = re.compile(r'[\s;]')
        if param.__class__.__name__ != "UniprotMapping":
            self.ownlog.msg(2, "Invalid parameter for read_mapping_uniprot()",
                            'ERROR')
            return {}
        mapping_o = {}
        mapping_i = {}
        scolend = re.compile(r'$;')
        rev = '' if not param.swissprot \
            else ' AND reviewed:%s' % param.swissprot
        query = 'organism:%u%s' % (int(ncbi_tax_id), rev)
        self.url = urls.urls['uniprot_basic']['url']
        self.post = {
            'query': query,
            'format': 'tab',
            'columns': 'id,%s%s' % (param.field, '' if param.subfield is None
                                    else '(%s)' % param.subfield)
        }
        self.url = '%s?%s' % (self.url, urllib.urlencode(self.post))
        c = curl.Curl(self.url, silent=False)
        data = c.result
        self.data = data
        data = [[[xx] if param.field == 'protein names' else [
            xxx for xxx in resep.split(scolend.sub('', xx.strip()))
            if len(xxx) > 0
        ] for xx in x.split('\t') if len(xx.strip()) > 0]
                for x in data.split('\n') if len(x.strip()) > 0]
        if len(data) > 0:
            del data[0]
            for l in data:
                if len(l) > 1:
                    l[1] = self.process_protein_name(l[1][0]) \
                        if param.field == 'protein names' else l[1]
                    for other in l[1]:
                        if other not in mapping_o:
                            mapping_o[other] = []
                        mapping_o[other].append(l[0][0])
                        if param.bi:
                            if l[0][0] not in mapping_i:
                                mapping_i[l[0][0]] = []
                            mapping_i[l[0][0]].append(other)
        self.mapping['to'] = mapping_o
        if param.bi:
            self.mapping['from'] = mapping_i

    def read_mapping_pickle(self, param, ncbi_tax_id = None):
        ncbi_tax_id = self.get_tax_id(ncbi_tax_id)
        if param.__class__.__name__ != "PickleMapping":
            self.ownlog.msg(2, "Invalid parameter for read_mapping_pickle()",
                            'ERROR')
            return False
        mapping = pickle.load(open(param.pickleFile, "rb"))
        if mapping.__class__.__name__ == "MappingTable":
            self = mapping
            return True
        else:
            return False

    def save_mapping_pickle(self, param):
        if param.__class__.__name__ != "PickleMapping":
            self.ownlog.msg(2, "Invalid parameter for save_mapping_pickle()",
                            'ERROR')
            return 1
        pickle.dump(self, open(param.pickleFile + ".pickle", "wb"))
        return 0

    def read_mapping_mysql(self, param):
        if param.mysql is None:
            self.ownlog.msg(2, 'No MySQL parameters given.', 'ERROR')
            return {"o": {}, "i": {}}
        tax_filter = ("" if param.ncbi_tax_id is None else
                      "AND %s = %u" % (param.ncbi_tax_id, self.ncbi_tax_id))
        query = """
            SELECT %s AS one,%s AS two FROM %s 
            WHERE %s IS NOT NULL AND %s IS NOT NULL %s""" % (
            param.fieldOne, param.fieldTwo, param.tableName, param.fieldOne,
            param.fieldTwo, tax_filter)
        try:
            param.mysql.run_query(query)
        except _mysql.Error as e:
            self.ownlog.msg(2, "MySQL error: %s\nFAILED QUERY: %s" %
                            (e, query), 'ERROR')
            return {"o": {}, "i": {}}
        total = len(param.mysql.result) + 1
        prg = progress.Progress(
            total=total, name="Processing data", interval=42)
        mapping_o = {}
        mapping_i = {}
        for rr in param.mysql.result:
            if rr["one"] not in mapping_o:
                mapping_o[rr["one"]] = []
            if rr["two"] not in mapping_i:
                mapping_i[rr["two"]] = []
            mapping_o[rr["one"]].append(rr["two"])
            mapping_i[rr["two"]].append(rr["one"])
            prg.step()
        self.mapping["to"] = mapping_o
        self.cleanDict(self.mapping["to"])
        if param.bi:
            self.mapping["from"] = mapping_i
            self.cleanDict(self.mapping["from"])
        prg.terminate()

    def process_protein_name(self, name):
        rebr = re.compile(r'\(([^\)]{3,})\)')
        resq = re.compile(r'\[([^\]]{3,})\]')
        names = [name.split('(')[0]]
        names += rebr.findall(name)
        others = flatList([x.split(';') for x in resq.findall(name)])
        others = [x.split(':')[1] if ':' in x else x for x in others]
        others = [x.split('(')[1] if '(' in x else x for x in others]
        names += others
        return [x.strip() for x in names]

    def id_max_len(self):
        if self.maxlOne is None:
            self.maxlOne = max(
                len(i) for i in flatList(self.mapping["to"].values()))
        if self.maxlTwo is None:
            self.maxlTwo = max(
                len(i) for i in flatList(self.mapping["from"].values()))
        return {"one": self.maxlOne, "two": self.maxlTwo}
    
    def get_tax_id(self, ncbi_tax_id):
        return ncbi_tax_id if ncbi_tax_id is not None else self.param.ncbi_tax_id


class Mapper(object):
    
    def __init__(self,
                 ncbi_tax_id=9606,
                 mysql_conf=(None, 'mapping'),
                 log=None,
                 cache=True,
                 cachedir='cache'):
        
        self.reup = re.compile(
            r'[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}'
        )
        self.cache = cache
        self.cachedir = cachedir
        if self.cache and not os.path.exists(self.cachedir):
            os.mkdir(self.cachedir)
        self.unmapped = []
        self.ownlog = log
        self.default_ncbi_tax_id = ncbi_tax_id
        self.tables = {}
        self.tables[self.default_ncbi_tax_id] = {}
        self.uniprot_mapped = []
        if log.__class__.__name__ != 'logw':
            self.session = common.gen_session_id()
            self.ownlog = logn.logw(self.session, 'INFO')
        else:
            self.ownlog = log
        self.mysql_conf = mysql_conf
        self.mysql = None
        self.trace = []
        self.name_types = {
            'uniprot_id': 'UniProtKB-ID',
            'embl': 'EMBL-CDS',
            'embl_id': 'EMBL',
            'entrez': 'GeneID',
            'gi': 'GI',
            'refseqp': 'RefSeq',
            'refseqn': 'RefSeq_NT',
            'ensembl': 'Ensembl',
            'ensg': 'Ensembl Genome',
            'ensp': 'Ensembl_PRO',
            'enst': 'Ensembl_TRS',
            'hgnc': 'HGNC'
        }
        self.types_name = dict(
            zip(self.name_types.values(), self.name_types.keys()))

    def reload(self):
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)

    def init_mysql(self):
        self.mysql = mysql.MysqlRunner(self.mysql_conf, log=self.ownlog)

    def which_table(self, nameType, targetNameType,
                    load=True, ncbi_tax_id = None):
        '''
        Returns the table which is suitable to convert an ID of
        nameType to targetNameType. If no such table have been loaded
        yet, it attempts to load from UniProt. If all attempts failed
        returns `None`.
        '''
        
        tbl = None
        ncbi_tax_id = self.get_tax_id(ncbi_tax_id)
        tblName = (nameType, targetNameType)
        tblNameRev = (targetNameType, nameType)
        
        if ncbi_tax_id not in self.tables:
            self.tables[ncbi_tax_id] = {}
        
        tables = self.tables[ncbi_tax_id]
        
        if tblName in tables:
            tbl = tables[tblName].mapping['to']
        
        elif tblNameRev in tables and \
                len(tables[tblNameRev].mapping['from']) > 0:
            tbl = tables[tblNameRev].mapping['from']
        
        elif load:
            
            for form in ['mapListUniprot', 'mapListBasic', 'mapListMirbase']:
                
                frm = getattr(maps, form)
                
                if tblName in frm:
                    self.load_mappings(maplst={tblName: frm[tblName]},
                                       ncbi_tax_id=ncbi_tax_id)
                    tbl = self.which_table(
                        nameType, targetNameType, load=False,
                        ncbi_tax_id = ncbi_tax_id)
                    break
                
                if tblNameRev in frm:
                    frm[tblNameRev].bi = True
                    self.load_mappings(maplst={tblNameRev: frm[tblNameRev]},
                                       ncbi_tax_id=ncbi_tax_id)
                    tbl = self.which_table(
                        nameType, targetNameType, load=False,
                        ncbi_tax_id = ncbi_tax_id)
                    break
                
                if tbl is not None:
                    break
            
            if tbl is None:
                
                if nameType in self.name_types:
                    this_param = input_formats.UniprotListMapping(
                        nameType = nameType,
                        targetNameType = targetNameType,
                        ncbi_tax_id = ncbi_tax_id)
                    
                    tables[tblName] = MappingTable(
                        nameType,
                        targetNameType,
                        this_param.typ,
                        this_param.__class__.__name__\
                            .replace('Mapping', '').lower(),
                        this_param,
                        ncbi_tax_id,
                        log=self.ownlog,
                        cache=self.cache,
                        cachedir=self.cachedir
                    )
                
                tbl = self.which_table(
                        nameType, targetNameType, load=False,
                        ncbi_tax_id = ncbi_tax_id)
            
            if tbl is None:
                
                if nameType in self.name_types:
                    self.load_uniprot_mappings([nameType])
                    
                    tbl = self.which_table(
                            nameType, targetNameType, load=False,
                            ncbi_tax_id = ncbi_tax_id)
        return tbl
    
    def get_tax_id(self, ncbi_tax_id):
        return self.default_ncbi_tax_id if ncbi_tax_id is None else ncbi_tax_id

    def map_name(self,
                 name,
                 nameType,
                 targetNameType,
                 ncbi_tax_id=None,
                 strict=False,
                 silent=True):
        r"""
        This function should be used to convert individual IDs.
        It takes care about everything, you don't need to think
        on the details. How does it work: looks up dictionaries 
        between the original and target ID type, if doesn't 
        find, attempts to load from the predefined inputs.
        If the original name is genesymbol, first it looks up
        among the preferred gene names from UniProt, if not 
        found, it takes an attempt with the alternative gene
        names. If the gene symbol still couldn't be found, and 
        strict = False, the last attempt only the first 5 chara-
        cters of the gene symbol matched. If the target name 
        type is uniprot, then it converts all the ACs to primary. 
        Then, for the Trembl IDs it looks up the preferred gene 
        names, and find Swissprot IDs with the same preferred 
        gene name.

        @name : str
            The original name which shall be converted.
        @nameType : str
            The type of the name.
            Available by default:
            - genesymbol (gene name)
            - entrez (Entrez Gene ID \[#\])
            - refseqp (NCBI RefSeq Protein ID \[NP\_\*|XP\_\*\])
            - ensp (Ensembl protein ID \[ENSP\*\])
            - enst (Ensembl transcript ID \[ENST\*\])
            - ensg (Ensembl genomic DNA ID \[ENSG\*\])
            - hgnc (HGNC ID \[HGNC:#\])
            - gi (GI number \[#\])
            - embl (DDBJ/EMBL/GeneBank CDS accession)
            - embl_id (DDBJ/EMBL/GeneBank accession)
            To use other IDs, you need to define the input method
            and load the table before calling :py:func:Mapper.map_name().

        """
        
        ncbi_tax_id = self.get_tax_id(ncbi_tax_id)
        if type(nameType) is list:
            mappedNames = []
            for nt in nameType:
                mappedNames += self.map_name(name, nt, targetNameType, strict,
                                             silent)
            return common.uniqList(mappedNames)
        if nameType == targetNameType:
            if targetNameType != 'uniprot':
                return [name]
            else:
                mappedNames = [name]
        elif nameType.startswith('refseq'):
            mappedNames = self.map_refseq(name,
                                          nameType,
                                          targetNameType,
                                          ncbi_tax_id = ncbi_tax_id,
                                          strict=strict)
        else:
            mappedNames = self._map_name(name,
                                         nameType,
                                         targetNameType,
                                         ncbi_tax_id)
        if not len(mappedNames):
            mappedNames = self._map_name(name.upper(),
                                         nameType,
                                         targetNameType,
                                         ncbi_tax_id)
        if not len(mappedNames) and \
            nameType not in set(['uniprot', 'trembl', 'uniprot-sec']):
            mappedNames = self._map_name(name.lower(),
                                         nameType,
                                         targetNameType,
                                         ncbi_tax_id)
        if not len(mappedNames) and nameType == 'genesymbol':
            mappedNames = self._map_name(name,
                                         'genesymbol-syn',
                                         targetNameType,
                                         ncbi_tax_id)
            if not strict and not len(mappedNames):
                mappedNames = self._map_name('%s1' % name,
                                             'genesymbol',
                                             targetNameType,
                                             ncbi_tax_id)
                if not len(mappedNames):
                    mappedNames = self._map_name(name,
                                                 'genesymbol5',
                                                 targetNameType,
                                                 ncbi_tax_id)
        
        if not len(mappedNames) and nameType == 'mir-mat-name':
            
            mappedNames = self._map_name(name,
                                         'mir-name',
                                         targetNameType,
                                         ncbi_tax_id)
        
        if targetNameType == 'uniprot':
            orig = mappedNames
            mappedNames = self.primary_uniprot(mappedNames)
            mappedNames = self.trembl_swissprot(mappedNames, ncbi_tax_id)
            if len(set(orig) - set(mappedNames)) > 0:
                self.uniprot_mapped.append((orig, mappedNames))
            mappedNames = [u for u in mappedNames if self.reup.match(u)]
        return common.uniqList(mappedNames)
    
    def map_names(self,
                 names,
                 nameType,
                 targetNameType,
                 ncbi_tax_id=None,
                 strict=False,
                 silent=True):
        """
        Same as `map_name` just with multiple IDs.
        
        """
        
        return (
            common.uniqList(
                itertools.chain(
                    *map(
                        lambda n:
                            self.map_name(n, nameType, targetNameType,
                                          ncbi_tax_id = ncbi_tax_id,
                                          strict = strict,
                                          silent = silent),
                        names
                    )
                )
            )
        )
    
    def map_refseq(self, refseq, nameType, targetNameType,
                   ncbi_tax_id, strict=False):
        mappedNames = []
        if '.' in refseq:
            mappedNames += self._map_name(refseq,
                                          nameType,
                                          targetNameType,
                                          ncbi_tax_id)
            if not len(mappedNames) and not strict:
                mappedNames += self._map_name(refseq.split('.')[0],
                                              nameType,
                                              targetNameType,
                                              ncbi_tax_id)
        if not len(mappedNames) and not strict:
            rstem = refseq.split('.')[0]
            for n in xrange(49):
                mappedNames += self._map_name('%s.%u' % (rstem, n),
                                              nameType,
                                              targetNameType,
                                              ncbi_tax_id)
        return mappedNames

    def _map_name(self, name, nameType, targetNameType, ncbi_tax_id):
        '''
        Once we have defined the name type and the target name type,
        this function looks it up in the most suitable dictionary.
        '''
        nameTypes = (nameType, targetNameType)
        nameTypRe = (targetNameType, nameType)
        tbl = self.which_table(nameType, targetNameType,
                               ncbi_tax_id = ncbi_tax_id)
        if tbl is None or name not in tbl:
            result = []
        elif name in tbl:
            result = tbl[name]
        # self.trace.append({'name': name, 'from': nameType, 'to': targetNameType,
        #    'result': result})
        return result

    def primary_uniprot(self, lst):
        '''
        For a list of UniProt IDs returns the list of primary ids.
        '''
        pri = []
        for u in lst:
            pr = self.map_name(u, 'uniprot-sec', 'uniprot-pri', ncbi_tax_id = 0)
            if len(pr) > 0:
                pri += pr
            else:
                pri.append(u)
        return list(set(pri))

    def trembl_swissprot(self, lst, ncbi_tax_id = None):
        '''
        For a list of Trembl and Swissprot IDs, returns possibly
        only Swissprot, mapping from Trembl to gene names, and 
        then back to Swissprot.
        '''
        ncbi_tax_id = self.get_tax_id(ncbi_tax_id)
        sws = []
        for tr in lst:
            sw = []
            gn = self.map_name(tr, 'trembl', 'genesymbol',
                               ncbi_tax_id = ncbi_tax_id)
            for g in gn:
                sw = self.map_name(g, 'genesymbol', 'swissprot',
                                   ncbi_tax_id = ncbi_tax_id)
            if len(sw) == 0:
                sws.append(tr)
            else:
                sws += sw
        sws = list(set(sws))
        return sws

    def has_mapping_table(self, nameTypeA, nameTypeB, ncbi_tax_id = None):
        ncbi_tax_id = self.get_tax_id(ncbi_tax_id)
        if (nameTypeA, nameTypeB) not in self.tables[ncbi_tax_id] and \
            ((nameTypeB, nameTypeA) not in self.tables[ncbi_tax_id] or
                len(self.tables[ncbi_tax_id][(nameTypeB, nameTypeA)].mapping['form']) == 0):
            self.map_table_error(nameTypeA, nameTypeB, ncbi_tax_id)
            return False
        else:
            return True

    def map_table_error(self, a, b, ncbi_tax_id):
        msg = ("Missing mapping table: from `%s` to `%s` at organism `%u` "\
            "mapping needed." % (a, b, ncbi_tax_id))
        sys.stdout.write(''.join(['\tERROR: ', msg, '\n']))
        self.ownlog.msg(2, msg, 'ERROR')

    def load_mappings(self, maplst=None, ncbi_tax_id = None):
        """
        mapList is a list of mappings to load;
        elements of mapList are dicts containing the 
        id names, molecule type, and preferred source
        e.g. ("one": "uniprot", "two": "refseq", "typ": "protein", 
        "src": "mysql", "par": "mysql_param/file_param")
        by default those are loaded from pickle files
        """
        
        ncbi_tax_id = self.get_tax_id(ncbi_tax_id)
        
        if maplst is None:
            try:
                maplst = maps.mapList
            except:
                self.ownlog.msg(1, 'load_mappings(): No input defined',
                                'ERROR')
                return None
        self.ownlog.msg(1, "Loading mapping tables...")
        
        for mapName, param in iteritems(maplst):
            
            tables = self.tables[ncbi_tax_id]
            
            param = param.set_organism(ncbi_tax_id)
            
            self.ownlog.msg(2, "Loading table %s ..." % str(mapName))
            sys.stdout.write("\t:: Loading '%s' to '%s' mapping table\n" %
                             (mapName[0], mapName[1]))
            
            typ = param.__class__.__name__
            
            if typ == 'FileMapping' and \
                    not os.path.isfile(param.input) and \
                    not hasattr(mapping_input, param.input):
                self.ownlog.msg(2, "Error: no such file: %s" % param.input,
                                "ERROR")
                continue
            
            if typ == 'PickleMapping' and \
                    not os.path.isfile(param.pickleFile):
                self.ownlog.msg(2, "Error: no such file: %s" %
                                m["par"].pickleFile, "ERROR")
                continue
            
            if typ == 'MysqlMapping':
                if not self.mysql:
                    self.ownlog.msg(2, "Error: no mysql server known.",
                                    "ERROR")
                    continue
                else:
                    if self.mysql is None:
                        self.init_mysql()
                    param.mysql = self.mysql
            
            if param.ncbi_tax_id != ncbi_tax_id:
                if typ == 'FileMapping':
                    sys.stdout.write('\t:: No translation table for organism `%u`. '\
                                     'Available for `%u` in file `%s`.\n' % \
                                     (ncbi_tax_id, param.ncbi_tax_id, param.input))
                    return None
            
            tables[mapName] = \
                MappingTable(
                    mapName[0],
                    mapName[1],
                    param.typ,
                    typ.replace('Mapping', '').lower(),
                    param,
                    ncbi_tax_id,
                    mysql=self.mysql,
                    log=self.ownlog,
                    cache=self.cache,
                    cachedir=self.cachedir
                )
            
            if ('genesymbol', 'uniprot') in tables \
                and ('genesymbol-syn', 'swissprot') in tables \
                and ('genesymbol5', 'uniprot') not in tables:
                self.genesymbol5(param.ncbi_tax_id)
            self.ownlog.msg(2, "Table %s loaded from %s." %
                            (str(mapName), param.__class__.__name__))

    def swissprots(self, lst):
        swprots = {}
        for u in lst:
            swprots[u] = self.map_name(u, 'uniprot', 'uniprot')
        return swprots

    def genesymbol5(self, ncbi_tax_id = None):
        ncbi_tax_id = self.get_tax_id(ncbi_tax_id)
        tables = self.tables[ncbi_tax_id]
        tables[('genesymbol5', 'uniprot')] = MappingTable(
            'genesymbol5',
            'uniprot',
            'protein',
            'genesymbol5',
            None,
            ncbi_tax_id,
            None,
            log=self.ownlog)
        tbl_gs5 = tables[('genesymbol5', 'uniprot')].mapping['to']
        tbls = [
            tables[('genesymbol', 'uniprot')].mapping['to'],
            tables[('genesymbol-syn', 'swissprot')].mapping['to']
        ]
        for tbl in tbls:
            for gs, u in iteritems(tbl):
                if len(gs) >= 5:
                    gs5 = gs[:5]
                    if gs5 not in tbl_gs5:
                        tbl_gs5[gs5] = []
                    tbl_gs5[gs5] += u
        tables[('genesymbol5', 'uniprot')].mapping['to'] = tbl_gs5

    def load_uniprot_mappings(self, ac_types=None, bi=False,
                              ncbi_tax_id = None):
        
        ncbi_tax_id = self.get_tax_id(ncbi_tax_id)
        tables = self.tables[ncbi_tax_id]
        ac_types = ac_types if ac_types is not None else self.name_types.keys()
        # creating empty MappingTable objects:
        for ac_typ in ac_types:
            tables[(ac_typ, 'uniprot')] = MappingTable(
                ac_typ,
                'uniprot',
                'protein',
                ac_typ,
                None,
                ncbi_tax_id,
                None,
                log=self.ownlog)
        # attempting to load them from Pickle
        i = 0
        for ac_typ in ac_types:
            md5ac = common.md5((ac_typ, 'uniprot', bi, ncbi_tax_id))
            cachefile = os.path.join('cache', md5ac)
            if self.cache and os.path.isfile(cachefile):
                tables[(ac_typ, 'uniprot')].mapping = \
                    pickle.load(open(cachefile, 'rb'))
                ac_types.remove(ac_typ)
                tables[(ac_typ, 'uniprot')].mid = md5ac
        # loading the remaining from the big UniProt mapping file:
        if len(ac_types) > 0:
            url = urls.urls['uniprot_idmap_ftp']['url']
            c = curl.Curl(url, silent=False, large=True)

            prg = progress.Progress(c.size, "Processing ID conversion list",
                                    99)
            for l in c.result:
                prg.step(len(l))
                l = l.decode('ascii').strip().split('\t')
                for ac_typ in ac_types:
                    if len(l) > 2 and self.name_types[ac_typ] == l[1]:
                        other = l[2].split('.')[0]
                        if l[2] not in tables[(ac_typ, 'uniprot'
                                                    )].mapping['to']:
                            tables[(
                                ac_typ, 'uniprot')].mapping['to'][other] = []
                        tables[(ac_typ, 'uniprot')].mapping['to'][other].\
                            append(l[0].split('-')[0])
                        if bi:
                            uniprot = l[0].split('-')[0]
                            if uniprot not in tables[(ac_typ, 'uniprot')].\
                                    mapping['from']:
                                tables[(ac_typ, 'uniprot')].\
                                    mapping['from'][uniprot] = []
                            tables[(ac_typ, 'uniprot')].mapping['from'][uniprot].\
                                append(other)
            prg.terminate()
            if self.cache:
                for ac_typ in ac_types:
                    md5ac = common.md5((ac_typ, bi))
                    cachefile = os.path.join('cache', md5ac)
                    pickle.dump(tables[(ac_typ, 'uniprot')].mapping,
                                open(cachefile, 'wb'))
    
    def save_all_mappings(self):
        self.ownlog.msg(1, "Saving all mapping tables...")
        for ncbi_tax_id in self.tables:
            for table in self.tables[ncbi_tax_id]:
                self.ownlog.msg(2, "Saving table %s ..." % table[0])
                param = mapping.pickleMapping(table)
                self.tables[m].save_mapping_pickle(param)
                self.ownlog.msg(2, "Table %s has been written to %s.pickle." %
                                (m, param.pickleFile))

    def load_uniprot_mapping(self, filename, ncbi_tax_id = None):
        '''
        This is a wrapper to load a ... mapping table.
        '''
        
        ncbi_tax_id = self.get_tax_id(ncbi_tax_id)
        tables = self.tables[ncbi_tax_id]
        umap = self.read_mapping_uniprot(filename, ncbi_tax_id,
                                         self.ownlog)
        for key, value in iteritems(umap):
            tables[key] = value

    def read_mapping_uniprot_mysql(self, filename, ncbi_tax_id, log, bi=False):
        if not os.path.isfile(filename):
            self.ownlog.msg(2, "No such file %s in read_mapping_uniprot()" %
                            (param, 'ERROR'))
        infile = codecs.open(filename, encoding='utf-8', mode='r')
        umap = {}
        self.ownlog.msg(2, "Loading UniProt mapping table from file %s" %
                        filename)
        for line in infile:
            if len(line) == 0:
                continue
            line = line.split()
            one = line[1].lower().replace("_", "-")
            uniprot = line[0]
            other = line[2]
            mapTableName = one + "_uniprot"
            if mapTableName not in umap:
                umap[mapTableName] = mappingTable(one, "uniprot", "protein",
                                                  "uniprot", None, None,
                                                  ncbi_tax_id, log)
            if other not in umap[mapTableName].mapping["to"]:
                umap[mapTableName].mapping["to"][other] = []
            umap[mapTableName].mapping["to"][other].append(uniprot)
            if bi:
                if uniprot not in umap[mapTableName].mapping["from"]:
                    umap[mapTableName].mapping["from"][uniprot] = []
                umap[mapTableName].mapping["from"][uniprot].append(other)
        self.ownlog.msg(2, "%u mapping tables from UniProt has been loaded" %
                        len(umap))
        return umap
