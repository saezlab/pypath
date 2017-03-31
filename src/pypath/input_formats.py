#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
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
import sys
import copy

__all__ = [
    'MysqlMapping', 'FileMapping', 'PickleMapping', 'ReadSettings', 'ReadList',
    'Reference', 'UniprotListMapping'
]


class MysqlMapping(object):
    def __init__(self,
                 tableName,
                 fieldOne,
                 fieldTwo,
                 db=None,
                 ncbi_tax_id=None,
                 bi=False,
                 mysql=None,
                 typ='protein'):
        self.tableName = tableName
        self.fieldOne = fieldOne
        self.fieldTwo = fieldTwo
        self.ncbi_tax_id = ncbi_tax_id
        self.db = db
        self.bi = bi
        self.mysql = mysql
        self.typ = typ
    
    def set_organism(self, ncbi_tax_id):
        other_organism = copy.deepcopy(self)
        other_organism.ncbi_tax_id = ncbi_tax_id
        return other_organism


class FileMapping(object):
    def __init__(self,
                 input,
                 oneCol,
                 twoCol,
                 separator=None,
                 header=0,
                 bi=False,
                 ncbi_tax_id=9606,
                 typ='protein'):
        
        self.input = input
        self.oneCol = oneCol
        self.twoCol = twoCol
        self.separator = separator
        self.header = header
        self.typ = typ
        self.bi = bi
        self.ncbi_tax_id = ncbi_tax_id
        self.inputArgs = {'organism': ncbi_tax_id}
    
    def set_organism(self, ncbi_tax_id):
        other_organism = copy.deepcopy(self)
        other_organism.ncbi_tax_id = ncbi_tax_id
        
        if 'organism' in other_organism.inputArgs:
            other_organism.inputArgs['organism'] = ncbi_tax_id
        
        return other_organism

# class UniprotMapping(object):

# def __init__(self, ac_type, field = None,
# bi = False, tax = 9606, swissprot = 'yes', subfield = None):
#self.bi = bi
#self.tax = int(tax)
#self.ac_type = ac_type
#self.typ = 'protein'
# if field is not None:
#self.field = field
#self.subfield = subfield
# elif self.ac_type in ac_types:
#self.field = ac_types[self.ac_type][0]
#self.subfield = ac_types[self.ac_type][1]
#self.swissprot = swissprot


class UniprotMapping(object):
    
    def __init__(self, nameType, bi=False, ncbi_tax_id=9606, swissprot='yes'):
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
        self.ncbi_tax_id = int(ncbi_tax_id)
        self.typ = 'protein'
        self.swissprot = swissprot
        self.nameType = nameType
        self.targetNameType = 'uniprot'
        self.field = None if nameType not in ac_query \
            else ac_query[nameType][0]
        self.subfield = None if nameType not in ac_query \
            else ac_query[nameType][1]
    
    def set_organism(self, ncbi_tax_id):
        other_organism = copy.deepcopy(self)
        other_organism.ncbi_tax_id = ncbi_tax_id
        return other_organism

class UniprotListMapping(object):
    
    def __init__(self, nameType, ac_name = None,
                 targetNameType = None, target_ac_name = None,
                 bi = False,
                 ncbi_tax_id = 9606,
                 swissprot = True):
        """
        Provides parameters for downloading mapping table from UniProt
        `Upload Lists` webservice.
        """
        self.swissprot = swissprot
        self.ac_mapping = ac_mapping
        self.bi = bi
        self.ncbi_tax_id = ncbi_tax_id
        self.typ = 'protein'
        self.nameType, self.ac_name = self.get_ac_type(nameType, ac_name)
        self.targetNameType, self.target_ac_name = \
            self.get_ac_type(targetNameType, target_ac_name)
    
    def get_ac_type(self, nameType, ac_name):
        nameType = nameType if nameType is not None else \
            ac_name if ac_name is not None else 'uniprot'
        if nameType not in self.ac_mapping and ac_name is None:
            sys.stdout.write('\t:: Unknown ID type: `%s`.\n' % nameType)
        else:
            ac_name = self.ac_mapping[nameType] if ac_name is None else ac_name
        return nameType, ac_name
    
    def set_organism(self, ncbi_tax_id):
        other_organism = copy.deepcopy(self)
        other_organism.ncbi_tax_id = ncbi_tax_id
        return other_organism

class PickleMapping(object):
    def __init__(self, pickleFile):
        self.pickleFile = pickleFile


class ReadSettings:
    def __init__(self,
                 name="unknown",
                 separator=None,
                 nameColA=0,
                 nameColB=1,
                 nameTypeA="uniprot",
                 nameTypeB="uniprot",
                 typeA="protein",
                 typeB="protein",
                 isDirected=False,
                 sign=False,
                 inFile=None,
                 references=False,
                 extraEdgeAttrs={},
                 extraNodeAttrsA={},
                 extraNodeAttrsB={},
                 header=False,
                 taxonA=9606,
                 taxonB=9606,
                 ncbiTaxId=False,
                 interactionType='PPI',
                 positiveFilters=[],
                 negativeFilters=[],
                 inputArgs={},
                 must_have_references=True,
                 huge=False,
                 resource=None):
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
        self.must_have_references = must_have_references and bool(references)
        self.huge = huge
        self.resource = self.name if resource is None else resource


class ReadList:
    def __init__(self,
                 name="unknown",
                 separator=None,
                 nameCol=0,
                 nameType="uniprot",
                 typ="protein",
                 inFile=None,
                 extraAttrs={},
                 header=False):
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
