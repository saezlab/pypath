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

import os
import sys
import codecs
import re
import urllib
import urllib2
import hashlib
import json
try:
    import cPickle as pickle
except:
    import pickle

# from bioigraph:
import progress
import logn
from common import *
import data_formats
import mysql
import dataio

__all__ = ['MappingTable', 'Mapper']

###
### functions to read and use mapping tables from UniProt, file, mysql or pickle 
###

class MappingTable(object):
    '''
    To initialize ID conversion tables for the first time
    data is downloaded from UniProt and read to dictionaries.
    It takes a couple of seconds. Data is saved to pickle 
    dumps, this way after tables load much faster.
    '''
    
    def __init__(self, one, two, typ, source, param, ncbi_tax_id, 
        mysql = None, log = None, cache = True, cachedir = 'cache'):
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
            self.session = gen_session_id()
            self.ownlog = logn.logw(self.session, 'INFO')
        else:
            self.ownlog = log
        if param is not None:
            self.mid = hashlib.md5(str((one, two, self.param.bi))).hexdigest()
            md5param = hashlib.md5(json.dumps(self.param.__dict__)).hexdigest()
            self.cachefile = os.path.join(self.cachedir, md5param)
            if self.cache and os.path.isfile(self.cachefile):
                self.mapping = pickle.load(open(self.cachefile, 'rb'))
            elif len(self.mapping['to']) == 0 or (param.bi and len(self.mapping['from']) == 0):
                if os.path.exists(self.cachefile):
                    os.remove(self.cachefile)
                if source == "mysql":
                    self.read_mapping_mysql(param)
                elif source == "file":
                    self.read_mapping_file(param)
                elif source == "pickle":
                    self.read_mapping_pickle(param)
                elif source == "uniprot":
                    self.read_mapping_uniprot(param)
                if len(self.mapping['to']) != 0 and (not param.bi or len(self.mapping['from']) != 0):
                    pickle.dump(self.mapping, open(self.cachefile, 'wb'))
    
    def cleanDict(self,mapping):
        for key, value in mapping.iteritems():
            mapping[key] = uniqList(value)
        return mapping
    
    def read_mapping_file(self, param):
        if param.__class__.__name__ != "FileMapping":
            self.ownlog.msg(2, "Invalid parameter for read_mapping_file()", 'ERROR')
            return {}
        if not os.path.exists(param.input) and not hasattr(dataio, param.input):
            return {}
        if hasattr(dataio, param.input):
            toCall = getattr(dataio, param.input)
            infile = toCall()
            total = sum([sys.getsizeof(i) for i in infile])
        else:
            infile = codecs.open(param.input, encoding = 'utf-8', mode = 'r')
            total = os.path.getsize(param.input)
        prg = progress.Progress(
            total = total, name = "Reading from file", interval = 18)
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
        if type(infile) is file:
            infile.close()
        self.mapping["to"] = mapping_o
        self.cleanDict(self.mapping["to"])
        if param.bi:
            self.mapping["from"] = mapping_i
            self.cleanDict(self.mapping["from"])
        prg.terminate()
    
    def read_mapping_uniprot(self, param):
        '''
        Downloads ID mappings directly from UniProt.
        See the names of possible identifiers here:
        http://www.uniprot.org/help/programmatic_access
        
        @param : UniprotMapping instance
        '''
        resep = re.compile(r'[\s;]')
        if param.__class__.__name__ != "UniprotMapping":
            self.ownlog.msg(2, "Invalid parameter for read_mapping_uniprot()",
                'ERROR')
            return {}
        mapping_o = {}
        mapping_i = {}
        scolend = re.compile(r'$;')
        rev = '' if param.swissprot is None \
            else ' AND reviewed:%s' % param.swissprot
        query = 'organism:%u%s' % (int(param.tax), rev)
        self.url = data_formats.urls['uniprot_basic']['url']
        self.post = {
            'query': query, 
            'format': 'tab', 
            'columns': 'id,%s%s' % (param.field, 
                '' if param.subfield is None else '(%s)'%param.subfield)
        }
        self.url = '%s?%s' % (self.url, urllib.urlencode(self.post))
        data = dataio.curl(self.url, silent = False)
        self.data = data
        data = [[[xx] if param.field == 'protein names' else \
            [xxx for xxx in resep.split(scolend.sub('', xx.strip())) if len(xxx) > 0] \
            for xx in x.split('\t') if len(xx.strip()) > 0] 
            for x in data.split('\n') if len(x.strip()) > 0]
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
    
    def read_mapping_pickle(self,param):
        if param.__class__.__name__ != "PickleMapping":
            self.ownlog.msg(2,"Invalid parameter for read_mapping_pickle()", 'ERROR')
            return False
        mapping = pickle.load( open( param.pickleFile, "rb" ) )
        if mapping.__class__.__name__ == "MappingTable":
            self = mapping
            return True
        else:
            return False
    
    def save_mapping_pickle(self,param):
        if param.__class__.__name__ != "PickleMapping":
            self.ownlog.msg(2,"Invalid parameter for save_mapping_pickle()", 'ERROR')
            return 1
        pickle.dump(self, open( param.pickleFile+".pickle", "wb" ))
        return 0
    
    def read_mapping_mysql(self, param):
        if param.mysql is None:
            self.ownlog.msg(2,'No MySQL parameters given.', 'ERROR')
            return {"o": {}, "i": {}}
        tax_filter = ("" if param.tax is None else "AND %s = %u" 
            % (param.tax, self.ncbi_tax_id))
        query = """
            SELECT %s AS one,%s AS two FROM %s 
            WHERE %s IS NOT NULL AND %s IS NOT NULL %s""" % (
                param.fieldOne, param.fieldTwo, param.tableName, 
                param.fieldOne, param.fieldTwo, tax_filter)
        try:
            param.mysql.run_query(query)
        except _mysql.Error, e:
            self.ownlog.msg(2,"MySQL error: %s\nFAILED QUERY: %s" % (e,query), 'ERROR')
            return {"o": {}, "i": {}}
        total = len(param.mysql.result) + 1
        prg = progress.Progress(
            total=total,name="Processing data",interval=42)
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
        self.mapping["to"] =  mapping_o
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
            self.maxlOne = max(len(i) for i in flatList(self.mapping["to"].values()))
        if self.maxlTwo is None:
            self.maxlTwo = max(len(i) for i in flatList(self.mapping["from"].values()))
        return {"one": self.maxlOne, "two": self.maxlTwo}

class Mapper(object):
    
    def __init__(self, ncbi_tax_id = 9606, mysql_conf = (None, 'mapping'), 
        log = None, cache = True, cachedir = 'cache'):
        self.cache = cache
        self.cachedir = cachedir
        if self.cache and not os.path.exists(self.cachedir):
            os.mkdir(self.cachedir)
        self.unmapped = []
        self.ownlog = log
        self.ncbi_tax_id = ncbi_tax_id
        self.tables = {}
        self.uniprot_mapped = []
        if log.__class__.__name__ != 'logw':
            self.session = gen_session_id()
            self.ownlog = logn.logw(self.session,'INFO')
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
            'ensg': 'Ensembl',
            'ensp': 'Ensembl_PRO',
            'enst': 'Ensembl_TRS',
            'hgnc': 'HGNC'
        }
        self.types_name = dict(zip(self.name_types.values(), self.name_types.keys()))
    
    def init_mysql(self):
        self.mysql = mysql.MysqlRunner(self.mysql_conf, log = self.ownlog)
    
    def which_table(self, nameType, targetNameType, load = True):
        '''
        Returns the table which is suitable to convert an ID of 
        nameType to targetNameType. If no such table have been loaded
        yet, it attempts to load from UniProt.
        '''
        tbl = None
        tblName = (nameType, targetNameType)
        tblNameRev = (targetNameType, nameType)
        if tblName in self.tables:
            tbl = self.tables[tblName].mapping['to']
        elif tblNameRev in self.tables and \
            len(self.tables[tblNameRev].mapping['from']) > 0:
            tbl = self.tables[tblNameRev].mapping['from']
        elif load:
            for form in ['mapListUniprot', 'mapListBasic']:
                frm = getattr(data_formats, form)
                if tblName in frm:
                    self.load_mappings(maps = {tblName: frm[tblName]})
                    tbl = self.which_table(nameType, targetNameType, load = False)
                    break
                if tblNameRev in frm:
                    frm[tblNameRev].bi = True
                    self.load_mappings(maps = {tblNameRev: frm[tblNameRev]})
                    tbl = self.which_table(nameType, targetNameType, load = False)
                    break
                if tbl is not None:
                    break
            if tbl is None:
                if nameType in self.name_types:
                    self.load_uniprot_mappings([nameType])
        return tbl
    
    def map_name(self, name, nameType, targetNameType, 
        strict = False, silent = True):
        '''
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
                genesymbol (gene name)
                entrez (Entrez Gene ID [#])
                refseqp (NCBI RefSeq Protein ID [NP_*|XP_*])
                ensp (Ensembl protein ID [ENSP*])
                enst (Ensembl transcript ID [ENST*])
                ensg (Ensembl genomic DNA ID [ENSG*])
                hgnc (HGNC ID [HGNC:#])
                gi (GI number [#])
                embl (DDBJ/EMBL/GeneBank CDS accession)
                embl_id (DDBJ/EMBL/GeneBank accession)
            To use other IDs, you need to define the input method
            and load the table before calling Mapper.map_name().
        '''
        # print '\tmapping %s' % name
        if nameType == targetNameType:
            if targetNameType != 'uniprot':
                return [ name ]
            else:
                mappedNames = [ name ]
        elif nameType.startswith('refseq'):
            mappedNames = self.map_refseq(name, nameType, targetNameType, strict = strict)
        else:
            mappedNames = self._map_name(name, nameType, targetNameType)
        if len(mappedNames) == 0:
            mappedNames = self._map_name(name.upper(), nameType, targetNameType)
        if len(mappedNames) == 0:
            mappedNames = self._map_name(name.lower(), nameType, targetNameType)
        if len(mappedNames) == 0 and nameType == 'genesymbol':
            mappedNames = self._map_name(name, 'genesymbol-syn', targetNameType)
            if not strict and len(mappedNames) == 0:
                mappedNames = self._map_name(name, 'genesymbol5', targetNameType)
        if targetNameType == 'uniprot':
            orig = mappedNames
            mappedNames = self.primary_uniprot(mappedNames)
            mappedNames = self.trembl_swissprot(mappedNames)
            if len(set(orig) - set(mappedNames)) > 0:
                self.uniprot_mapped.append((orig, mappedNames))
        # print '\tmapped to %s' % str(mappedNames)
        return list(set(mappedNames))
    
    def map_refseq(self, refseq, nameType, targetNameType, strict = False):
        mappedNames = []
        if '.' in refseq:
            mappedNames += self._map_name(refseq, nameType, targetNameType)
            if len(mappedNames) == 0 and not strict:
                mappedNames += self._map_name(refseq.split('.')[0], nameType, targetNameType)
        if len(mappedNames) == 0 and not strict:
            rstem = refseq.split('.')[0]
            for n in xrange(49):
                mappedNames += self._map_name('%s.%u'%(rstem, n), nameType, targetNameType)
        return mappedNames
    
    def _map_name(self, name, nameType, targetNameType):
        '''
        Once we have defined the name type and the target name type,
        this function looks it up in the most suitable dictionary.
        '''
        nameTypes = (nameType, targetNameType)
        nameTypRe = (targetNameType, nameType)
        tbl = self.which_table(nameType, targetNameType)
        if tbl is None or name not in tbl:
            result = []
        elif name in tbl:
            result = tbl[name]
        #self.trace.append({'name': name, 'from': nameType, 'to': targetNameType, 
        #    'result': result})
        return result
    
    def primary_uniprot(self, lst):
        '''
        For a list of UniProt IDs returns the list of primary ids.
        '''
        pri = []
        for u in lst:
            pr = self.map_name(u, 'uniprot-sec', 'uniprot-pri')
            if len(pr) > 0:
                pri += pr
            else:
                pri.append(u)
        return list(set(pri))
    
    def trembl_swissprot(self, lst):
        '''
        For a list of Trembl and Swissprot IDs, returns possibly
        only Swissprot, mapping from Trembl to gene names, and 
        then back to Swissprot.
        '''
        sws = []
        for tr in lst:
            sw = []
            gn = self.map_name(tr, 'trembl', 'genesymbol')
            for g in gn:
                sw = self.map_name(g, 'genesymbol', 'swissprot')
            if len(sw) == 0:
                sws.append(tr)
            else:
                sws += sw
        sws = list(set(sws))
        return sws
    
    def has_mapping_table(self, nameTypeA, nameTypeB):
        if (nameTypeA, nameTypeB) not in self.tables and \
            ((nameTypeB, nameTypeA) not in self.tables or \
                len(self.tables[(nameTypeB, nameTypeA)].mapping['form']) == 0):
            self.map_table_error(nameTypeA, nameTypeB)
            return False
        else:
            return True
    
    def map_table_error(self, a, b):
        msg = ("Missing mapping table: from %s to %s mapping needed." % (a, b))
        sys.stdout.write(''.join(['\tERROR: ',msg,'\n']))
        self.ownlog.msg(2,msg,'ERROR')
    
    def load_mappings(self, maps = None):
        '''
        mapList is a list of mappings to load;
        elements of mapList are dicts containing the 
        id names, molecule type, and preferred source
        e.g. ("one": "uniprot", "two": "refseq", "typ": "protein", 
        "src": "mysql", "par": "mysql_param/file_param")
        by default those are loaded from pickle files
        '''
        if maps is None: 
            try:
                maps = data_formats.mapList
            except:
                self.ownlog.msg(1, 'load_mappings(): No input defined','ERROR')
                return None
        self.ownlog.msg(1, "Loading mapping tables...")
        for mapName, param in maps.iteritems():
            self.ownlog.msg(2, "Loading table %s ..." % str(mapName))
            sys.stdout.write("\t:: Loading '%s' to '%s' mapping table\n" % \
                (mapName[0], mapName[1]))
            if param.__class__.__name__ == 'FileMapping' and \
                not os.path.isfile(param.input) and \
                not hasattr(dataio, param.input):
                self.ownlog.msg(2,"Error: no such file: %s" % \
                    param.input, "ERROR")
                continue
            if param.__class__.__name__ == 'PickleMapping' and \
                not os.path.isfile(param.pickleFile):
                self.ownlog.msg(2,"Error: no such file: %s" % \
                    m["par"].pickleFile, "ERROR")
                continue
            if param.__class__.__name__ == 'MysqlMapping':
                if not self.mysql:
                    self.ownlog.msg(2,"Error: no mysql server known.", "ERROR")
                    continue
                else:
                    if self.mysql is None:
                        self.init_mysql()
                    param.mysql = self.mysql
            self.tables[mapName] = MappingTable(
                mapName[0],
                mapName[1],
                param.typ,
                param.__class__.__name__.replace('Mapping', '').lower(),
                param,
                self.ncbi_tax_id,
                mysql = self.mysql,
                log = self.ownlog,
                cache = self.cache,
                cachedir = self.cachedir)
            if ('genesymbol', 'uniprot') in self.tables and \
                ('genesymbol-syn', 'swissprot') in self.tables and \
                ('genesymbol5', 'uniprot') not in self.tables:
                self.genesymbol5()
            self.ownlog.msg(2, "Table %s loaded from %s." % (str(mapName), \
                param.__class__.__name__))
    
    def swissprots(self,lst):
        swprots = {}
        for u in lst:
            swprots[u] = self.map_name(u, 'uniprot', 'uniprot')
        return swprots
    
    def genesymbol5(self):
        self.tables[('genesymbol5', 'uniprot')] = MappingTable(
            'genesymbol5', 
            'uniprot',
            'protein',
            'genesymbol5',
            None,
            self.ncbi_tax_id, 
            None,
            log = self.ownlog)
        tbl_gs5 = self.tables[('genesymbol5', 'uniprot')].mapping['to']
        tbls = [self.tables[('genesymbol', 'uniprot')].mapping['to'], 
            self.tables[('genesymbol-syn', 'swissprot')].mapping['to']]
        for tbl in tbls:
            for gs, u in tbl.iteritems():
                if len(gs) >= 5:
                    gs5 = gs[:5]
                    if gs5 not in tbl_gs5:
                        tbl_gs5[gs5] = []
                    tbl_gs5[gs5] += u
        self.tables[('genesymbol5', 'uniprot')].mapping['to'] = tbl_gs5
    
    def load_uniprot_mappings(self, ac_types = None, bi = False):
        ac_types = ac_types if ac_types is not None else self.name_types.keys()
        # creating empty MappingTable objects:
        for ac_typ in ac_types:
            self.tables[(ac_typ, 'uniprot')] = MappingTable(
                ac_typ, 
                'uniprot',
                'protein',
                ac_typ,
                None, 
                self.ncbi_tax_id, 
                None,
                log = self.ownlog)
        # attempting to load them from Pickle
        i = 0
        for ac_typ in ac_types:
            md5ac = hashlib.md5(str((ac_typ, 'uniprot', bi))).hexdigest()
            cachefile = os.path.join('cache', md5ac)
            if self.cache and os.path.isfile(cachefile):
                self.tables[(ac_typ, 'uniprot')].mapping = \
                    pickle.load(open(cachefile, 'rb'))
                ac_types.remove(ac_typ)
                self.tables[(ac_typ, 'uniprot')].mid = md5ac
        # loading the remaining from the big UniProt mapping file:
        if len(ac_types) > 0:
            url = data_formats.urls['uniprot_idmap_ftp']['url']
            data = dataio.curl(url, silent = False)
            data = data.split('\n')
            prg = progress.Progress(len(data), "Processing ID conversion list", 99)
            for l in data:
                prg.step()
                l = l.split('\t')
                for ac_typ in ac_types:
                    if len(l) > 2 and self.name_types[ac_typ] == l[1]:
                        other = l[2].split('.')[0]
                        if l[2] not in self.tables[(ac_typ, 'uniprot')].mapping['to']:
                            self.tables[(ac_typ, 'uniprot')].mapping['to'][other] = []
                        self.tables[(ac_typ, 'uniprot')].mapping['to'][other].\
                            append(l[0].split('-')[0])
                        if bi:
                            uniprot = l[0].split('-')[0]
                            if uniprot not in self.tables[(ac_typ, 'uniprot')].\
                                mapping['from']:
                                self.tables[(ac_typ, 'uniprot')].\
                                mapping['from'][uniprot] = []
                            self.tables[(ac_typ, 'uniprot')].mapping['from'][uniprot].\
                                append(other)
            prg.terminate()
            if self.cache:
                for ac_typ in ac_types:
                    md5ac = hashlib.md5(str((ac_typ, bi))).hexdigest()
                    cachefile = os.path.join('cache', md5ac)
                    pickle.dump(self.tables[(ac_typ, 'uniprot')].mapping, 
                        open(cachefile, 'wb'))
    
    def save_all_mappings(self):
        self.ownlog.msg(1, "Saving all mapping tables...")
        for m in self.tables:
            self.ownlog.msg(2, "Saving table %s ..." % m)
            param = mapping.pickleMapping(m)
            self.tables[m].save_mapping_pickle(param)
            self.ownlog.msg(2, "Table %s has been written to %s.pickle." % (
                m, param.pickleFile))
    
    def load_uniprot_mapping(self,filename):
        '''
        This is a wrapper to load a ... mapping table.
        '''
        umap = self.read_mapping_uniprot(filename, self.ncbi_tax_id, self.ownlog)
        for key, value in umap.iteritems():
            self.tables[key] = value

    def read_mapping_uniprot_mysql(self, filename, ncbi_tax_id, log, bi = False):
        if not os.path.isfile(filename):
            self.ownlog.msg(2,"No such file %s in read_mapping_uniprot()" % (
                param, 'ERROR'))
        infile = codecs.open(filename, encoding='utf-8', mode='r')
        umap = {}
        self.ownlog.msg(2,"Loading UniProt mapping table from file %s" % filename)
        for line in infile:
            if len(line) == 0:
                continue
            line = line.split()
            one = line[1].lower().replace("_", "-")
            uniprot = line[0]
            other = line[2]
            mapTableName = one + "_uniprot"
            if mapTableName not in umap:
                umap[mapTableName] = mappingTable(
                    one, "uniprot",
                    "protein", "uniprot",
                    None,None, ncbi_tax_id, log)
            if other not in umap[mapTableName].mapping["to"]:
                umap[mapTableName].mapping["to"][other] = []
            umap[mapTableName].mapping["to"][other].append(uniprot)
            if bi:
                if uniprot not in umap[mapTableName].mapping["from"]:
                    umap[mapTableName].mapping["from"][uniprot] = []
                umap[mapTableName].mapping["from"][uniprot].append(other)
        self.ownlog.msg(2,"%u mapping tables from UniProt has been loaded" % len(umap))
        return umap
