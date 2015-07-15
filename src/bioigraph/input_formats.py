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

import codecs

__all__ = ['MysqlMapping','FileMapping','ReferenceList',
           'PickleMapping','ReadSettings','ReadList']

class MysqlMapping(object):
        
    def __init__(self, tableName, fieldOne, fieldTwo, db = None, tax = None, 
            bi = False, mysql = None, typ = 'protein'):
        self.tableName = tableName
        self.fieldOne = fieldOne
        self.fieldTwo = fieldTwo
        self.tax = tax
        self.db = db
        self.bi = bi
        self.mysql = mysql
        self.typ = typ

class FileMapping(object):
        
    def __init__(self, input, oneCol, twoCol, separator = None,
        header = 0, bi = False, tax = 9606, typ = 'protein'):
        self.input = input
        self.oneCol = oneCol
        self.twoCol = twoCol
        self.separator = separator
        self.header = header
        self.typ = typ
        self.bi = bi

#class UniprotMapping(object):
    
    #def __init__(self, ac_type, field = None, 
        #bi = False, tax = 9606, swissprot = 'yes', subfield = None):
        #self.bi = bi
        #self.tax = int(tax)
        #self.ac_type = ac_type
        #self.typ = 'protein'
        #if field is not None:
            #self.field = field
            #self.subfield = subfield
        #elif self.ac_type in ac_types:
            #self.field = ac_types[self.ac_type][0]
            #self.subfield = ac_types[self.ac_type][1]
        #self.swissprot = swissprot

class UniprotMapping(object):
    
    def __init__(self, nameType, bi = False, tax = 9606, swissprot = 'yes'):
        '''
        Defines an ID conversion table to retrieve from UniProt.
        
        @nameType : str
            Type of accession numbers you would like to translate.
        @targetNameType : str
            Type of accession numbers you would like to translate to.
        @bi : bool
            Build the mapping table only from original AC to target AC,
            or if bi = True, the reverse table is also generated (from 
            target to original). 
        @tax : int
            NCBI Taxonomy ID of the organism of interest.
        @swissprot : str
            Look for SwissProt or Trembl.
            Passed directly to UniProt`s `reviewed` parameter. `yes` or `no`
            To fetch Trembl and SwissProt together, set value to None.
        @mapping : bool
            Get the data from UniProt`s programmatic access query interface,
            (uniprot.org/uniprot) or the batch retrieval/id mapping service 
            (uniprot.org/mapping). These have slightly different APIs and 
            capabilities. Some IDs can be obtained from the former, some 
            from the latter.
        '''
        self.bi = bi
        self.tax = int(tax)
        self.typ = 'protein'
        self.swissprot = swissprot
        self.nameType = nameType
        self.targetNameType = 'uniprot'
        self.field = None if nameType not in ac_query \
            else ac_query[nameType][0]
        self.subfield = None if nameType not in ac_query \
            else ac_query[nameType][1]

class ReferenceList(object):
    
    def __init__(self,nameType,typ,tax,inFile):
        self.infile = inFile
        self.nameType = nameType
        self.typ = typ
        self.tax = tax
    
    def load(self):
        f = codecs.open(self.infile,encoding='utf-8',mode='r')
        lst = []
        for l in f:
            lst.append(l.strip())
        f.close()
        self.lst = set(lst)

class PickleMapping(object):
    
    def __init__(self, pickleFile):
        self.pickleFile = pickleFile

class ReadSettings:
    
    def __init__(self, name = "unknown", separator = None, nameColA = 0, nameColB = 1, 
            nameTypeA = "uniprot", nameTypeB = "uniprot", typeA = "protein", 
            typeB = "protein", isDirected = False, sign = False, inFile = None, 
            references = False, extraEdgeAttrs = {}, extraNodeAttrsA = {}, 
            extraNodeAttrsB = {}, header = False, taxonA = False, taxonB = False, 
            ncbiTaxId = False, interactionType = 'PPI', 
            positiveFilters = [], negativeFilters = [], inputArgs = {}):
        self.typeA = typeA
        self.typeB = typeB
        self.nameColA = nameColA
        self.nameColB = nameColB
        self.nameTypeA = nameTypeA
        self.nameTypeB = nameTypeB
        self.isDirected = isDirected
        self.inFile = inFile
        self.extraEdgeAttrs = extraEdgeAttrs
        self.extraNodeAttrsA = extraNodeAttrsA
        self.extraNodeAttrsB = extraNodeAttrsB
        self.name = name
        self.separator = separator
        self.header = header
        self.refs = references
        self.sign = sign
        self.taxonA = taxonA
        self.taxonB = taxonB
        self.ncbiTaxId = ncbiTaxId
        self.intType = interactionType
        self.positiveFilters = positiveFilters
        self.negativeFilters = negativeFilters
        self.inputArgs = inputArgs

class ReadList:
    
    def __init__(self, name="unknown", separator=None, nameCol=0,
            nameType = "uniprot", typ="protein", inFile=None,
            extraAttrs={},header=False):
        self.typ = typ
        self.nameCol = nameCol
        self.nameType = nameType
        self.inFile = inFile
        self.extraAttrs = extraAttrs
        self.name = name
        self.separator = separator
        self.header = header

ac_query = {
    'genesymbol': ['genes', 'PREFERRED'],
    'genesymbol-syn': ['genes', 'ALTERNATIVE'],
    'hgnc': ['database', 'HGNC'],
    'embl': ['database', 'embl'],
    'entrez': ['database', 'geneid'],
    'refseqp': ['database', 'refseq'],
    'enst': ['database', 'ensembl'],
    'uniprot-entry': ['entry name', None],
    'protein-name': ['protein names', None]
}

ac_mapping = {
    'uniprot': 'ACC',
    'uniprot_id': 'ID',
    'embl': 'EMBL',
    'embl_id': 'EMBL_ID',
    'pir': 'PIR',
    'entrez': 'P_ENTREZGENEID',
    'gi': 'P_GI',
    'refseqp': 'P_REFSEQ_AC',
    'refseqn': 'REFSEQ_NT_ID',
    'ensembl': 'ENSEMBL_ID',
    'ensp': 'ENSEMBL_PRO_ID',
    'enst': 'ENSEMBL_TRS_ID',
    'ensg': 'ENSEMBLGENOME_ID',
    'ensgp': 'ENSEMBLGENOME_PRO_ID',
    'ensgt': 'ENSEMBLGENOME_TRS_ID',
    'hgnc': 'HGNC_ID'
}
