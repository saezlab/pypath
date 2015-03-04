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
import _mysql
import progress
import logn
from common import *
import data_formats
import mysql

__all__ = ['MappingTable','Mapper']

###
### functions to read and use mapping tables from textfile, mysql or pickle 
###

class MappingTable(object):
    
    def __init__(self,one,two,typ,source,param,ncbi_tax_id,mysql=None,log=None):
        self.one = one
        self.two = two
        self.typ = typ
        self.maxlOne = None
        self.maxlTwo = None
        self.ncbi_tax_id = ncbi_tax_id
        self.mysql = mysql
        self.mapping = {"to": {}, "from": {}}
        if log.__class__.__name__ != 'logw':
            self.session = gen_session_id()
            self.ownlog = logn.logw(self.session,'INFO')
        else:
            self.ownlog = log
        if source == "mysql":
            self.read_mapping_mysql(param)
        elif source == "file":
            self.read_mapping_file(param)
        elif source == "pickle":
            self.read_mapping_pickle(param)
    
    def cleanDict(self,mapping):
        for key, value in mapping.iteritems():
            mapping[key] = uniqList(value)
        return mapping
    
    def read_mapping_file(self,param):
        if param.__class__.__name__ != "FileMapping":
            self.ownlog.msg(2,"Invalid parameter for read_mapping_file()", 'ERROR')
            return {}
        total = os.path.getsize(param.textFile)
        prg = progress.Progress(
            total=total,name="Reading from file",interval=18)
        infile = codecs.open(param.textFile, encoding='utf-8', mode='r')
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
            prg.step(len(line))
            line = line.rstrip().split(param.separator)
            if len(line) > max([param.oneCol, param.twoCol]):
                if line[param.oneCol] not in mapping_o:
                    mapping_o[line[param.oneCol]] = []
                if line[param.twoCol] not in mapping_i:
                    mapping_i[line[param.twoCol]] = []
                mapping_o[line[param.oneCol]].append(line[param.twoCol])
                mapping_i[line[param.twoCol]].append(line[param.oneCol])
        infile.close()
        self.mapping["to"] =  mapping_o
        self.cleanDict(self.mapping["to"])
        if param.bi:
            self.mapping["from"] = mapping_i
            self.cleanDict(self.mapping["from"])
        prg.terminate()
    
    def read_mapping_pickle(self,param):
        if param.__class__.__name__ != "PickleMapping":
            self.ownlog.msg(2,"Invalid parameter for read_mapping_pickle()", 'ERROR')
            return False
        mapping = pickle.load( open( param.pickleFile, "rb" ) )
        if mapping.__class__.__name__ == "mappingTable":
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
    
    def read_mapping_mysql(self,param):
        if self.mysql is None:
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
            self.mysql.run_query(query)
        except _mysql.Error, e:
            self.ownlog.msg(2,"MySQL error: %s\nFAILED QUERY: %s" % (e,query), 'ERROR')
            return {"o": {}, "i": {}}
        total = len(self.mysql.result) + 1
        prg = progress.Progress(
            total=total,name="Processing data",interval=42)
        mapping_o = {}
        mapping_i = {}
        for rr in self.mysql.result:
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
    
    def id_max_len(self):
        if self.maxlOne is None or self.maxlTwo is None:
            maxlOne = 0
            maxlTwo = 0
            for i in self.mapping["to"].values():
                for j in i:
                    if len(j) > maxlOne:
                        maxlOne = len(j)
            for i in self.mapping["from"].values():
                for j in i:
                    if len(j) > maxlTwo:
                        maxlTwo = len(j)
            self.maxlOne = maxlOne
            self.maxlTwo = maxlTwo
        return {"one": self.maxlOne, "two": self.maxlTwo}

class Mapper(object):
    
    def __init__(self,ncbi_tax_id,mysql_conf=(None,'mapping'),log=None):
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
        self.mysql = mysql.MysqlRunner(self.mysql_conf,log=self.ownlog)
    
    def map_name(self, name, nameType, targetNameType):
        if nameType == targetNameType:
            if targetNameType != 'uniprot':
                return [ name ]
            else:
                mappedNames = [ name ]
        elif (nameType+'_'+targetNameType in self.tables and 
            name in self.tables[nameType+'_'+targetNameType].mapping['to']):
            mappedNames = self.tables[nameType+'_'+targetNameType].mapping['to'][name]
        elif (targetNameType+'_'+nameType in self.tables and 
            name in self.tables[targetNameType+'_'+nameType].mapping['from']):
            mappedNames = self.tables[targetNameType+'_'+nameType].mapping['from'][name]
        elif (nameType+'-fallback_'+targetNameType in self.tables and 
            name in self.tables[nameType+'-fallback_'+targetNameType].mapping['to']):
            mappedNames = self.tables[nameType + '-fallback_' + \
                targetNameType].mapping['to'][name]
        else:
            if name not in self.unmapped:
                self.unmapped.append(name)
            return [ 'unmapped' ]
        if targetNameType == 'uniprot':
            orig = mappedNames
            mappedNames = self.get_primary_uniprot(mappedNames)
            mappedNames = self.trembl_swissprot(mappedNames)
            if len(set(orig)-set(mappedNames)) > 0:
                self.uniprot_mapped.append((orig,mappedNames))
        return list(set(mappedNames))
    
    def get_primary_uniprot(self, lst):
        # for a list of UniProt IDs returns the list of primary ids
        pri = []
        for u in lst:
            if u in self.tables["uniprot-sec_uniprot-pri"].mapping["to"]:
                pri += self.tables["uniprot-sec_uniprot-pri"].mapping["to"][u]
            else:
                pri.append(u)
        return list(set(pri))
    
    def trembl_swissprot(self, lst):
        # for a list of Trembl and Swissprot IDs, returns possibly
        # only Swissprot, mapping through gene names
        if not self.has_mapping_table("trembl", "genesymbol"):
            return None
        if not self.has_mapping_table("genesymbol", "swissprot"):
            return None
        sw = []
        for t in lst:
            if t in self.tables["trembl_genesymbol"].mapping["to"]:
                gn = self.tables["trembl_genesymbol"].mapping["to"][t]
                for g in gn:
                    if g in self.tables["genesymbol_swissprot"].mapping["to"]:
                        sw += self.tables["genesymbol_swissprot"].mapping["to"][g]
            else:
                sw.append(t)
        sw = list(set(sw))
        return sw
    
    def has_mapping_table(self,frm,to):
        if '_'.join([frm,to]) not in self.tables:
            self.map_table_error(frm, to)
            return False
        else:
            return True
    
    def map_table_error(self,a,b):
        msg = ("Missing mapping table: from %s to %s mapping needed." % (a, b))
        sys.stdout.write(''.join(['\tERROR: ',msg,'\n']))
        self.ownlog.msg(2,msg,'ERROR')
    
    def load_mappings(self,maps=None):
        if maps is None:
            try:
                maps = data_formats.mapList
            except:
                self.ownlog.msg(1, 'load_mappings(): No input defined','ERROR')
                return None
        self.ownlog.msg(1, "Loading mapping tables...")
        # mapList is a list of desired mappings;
        # elements of mapList are dicts containing the 
        # id names, molecule type, and preferred source
        # e.g. ("one": "uniprot", "two": "refseq", "typ": "protein", 
        # "src": "mysql", "par": "mysql_param/file_param")
        # by default those are loaded from pickle files
        for m in maps:
            mapName = m["one"]+"_"+m["two"]
            self.ownlog.msg(2, "Loading table %s ..." % mapName)
            sys.stdout.write("Loading '%s' to '%s' mapping table\n" % (m['one'],m['two']))
            if m["src"] == "file" and not os.path.isfile(m["par"].textFile):
                self.ownlog.msg(2,"Error: no such file: %s" % m["par"].textFile, "ERROR")
                continue
            if m["src"] == "pickle" and not os.path.isfile(m["par"].pickleFile):
                self.ownlog.msg(2,"Error: no such file: %s" % m["par"].pickleFile, "ERROR")
                continue
            if m["src"] == "mysql" and not self.mysql:
                self.ownlog.msg(2,"Error: no mysql server known.", "ERROR")
                continue
            self.tables[mapName] = MappingTable(
                m["one"],
                m["two"],
                m["typ"],
                m["src"],
                m['par'],
                self.ncbi_tax_id,
                mysql = self.mysql,
                log = self.ownlog)
            self.ownlog.msg(2, "Table %s loaded from %s." % (mapName, m['src']))
    
    def swissprots(self,lst):
        swprots = {}
        for u in lst:
            swprots[u] = self.map_name(u,'uniprot','uniprot')
        return swprots
    
    def save_all_mappings(self):
        self.ownlog.msg(1, "Saving all mapping tables...")
        for m in self.tables:
            self.ownlog.msg(2, "Saving table %s ..." % m)
            param = mapping.pickleMapping(m)
            self.tables[m].save_mapping_pickle(param)
            self.ownlog.msg(2, "Table %s has been written to %s.pickle." % (
                m, param.pickleFile))
    
    def load_uniprot_mapping(self,filename):
        # this is a wrapper to load the mapping table
        # downloaded from uniprot directly
        umap = self.read_mapping_uniprot(filename, self.ncbi_tax_id, self.ownlog)
        for key, value in umap.iteritems():
            self.tables[key] = value

    def read_mapping_uniprot(self,filename, ncbi_tax_id, log, bi=False):
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

#def read_uniprot_primary(filename, log):
#    
