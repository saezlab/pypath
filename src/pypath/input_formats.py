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
#  Website: http://pypath.omnipathdb.org/
#

import codecs
import sys
import copy

__all__ = [
    'MysqlMapping', 'FileMapping', 'PickleMapping', 'ReadSettings', 'ReadList',
    'Reference', 'UniprotListMapping'
]


class FileMapping(object):
    
    def __init__(
            self,
            input,
            col_a,
            col_b,
            separator = None,
            header = 0,
            bi_directional  =  False,
            ncbi_tax_id = 9606,
            entity_type = 'protein',
        ):
        
        self.input = input
        self.col_a = col_a
        self.col_b = col_b
        self.separator = separator
        self.header = header
        self.entity_type = entity_type
        self.bi_directional = bi_directional
        self.ncbi_tax_id = ncbi_tax_id
        self.input_args = {'organism': ncbi_tax_id}
    
    
    def set_organism(self, ncbi_tax_id):
        
        other_organism = copy.deepcopy(self)
        other_organism.ncbi_tax_id = ncbi_tax_id
        
        if 'organism' in other_organism.input_args:
            
            other_organism.input_args['organism'] = ncbi_tax_id
        
        return other_organism


class UniprotMapping(object):

    def __init__(
            self,
            name_type,
            bi_directional = False,
            ncbi_tax_id = 9606,
            swissprot = 'yes',
        ):
        '''
        Defines an ID conversion table to retrieve from UniProt.

        @name_type : str
            Type of accession numbers you would like to translate.
        @target_name_type : str
            Type of accession numbers you would like to translate to.
        @bi_directional : bool
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
        self.bi_directional = bi_directional
        self.ncbi_tax_id = int(ncbi_tax_id)
        self.typ = 'protein'
        self.swissprot = swissprot
        self.name_type = name_type
        self.target_name_type = 'uniprot'
        self.field = None if name_type not in ac_query \
            else ac_query[name_type][0]
        self.subfield = None if name_type not in ac_query \
            else ac_query[name_type][1]
    
    
    def set_organism(self, ncbi_tax_id):
        
        other_organism = copy.deepcopy(self)
        other_organism.ncbi_tax_id = ncbi_tax_id
        return other_organism

class UniprotListMapping(object):

    def __init__(
            self,
            name_type,
            ac_name = None,
            target_name_type = None,
            target_ac_name = None,
            bi_directional = False,
            ncbi_tax_id = 9606,
            swissprot = True,
        ):
        """
        Provides parameters for downloading mapping table from UniProt
        `Upload Lists` webservice.
        """
        self.swissprot = swissprot
        self.ac_mapping = ac_mapping
        self.bi = bi
        self.ncbi_tax_id = ncbi_tax_id
        self.typ = 'protein'
        self.name_type, self.ac_name = self.get_ac_type(name_type, ac_name)
        self.target_name_type, self.target_ac_name = \
            self.get_ac_type(target_name_type, target_ac_name)
    
    
    def get_ac_type(self, name_type, ac_name):
        
        name_type = name_type if name_type is not None else \
            ac_name if ac_name is not None else 'uniprot'
        if name_type not in self.ac_mapping and ac_name is None:
            sys.stdout.write('\t:: Unknown ID type: `%s`.\n' % name_type)
        else:
            ac_name = (
                self.ac_mapping[name_type] if ac_name is None else ac_name
            )
        return name_type, ac_name
    
    
    def set_organism(self, ncbi_tax_id):
        
        other_organism = copy.deepcopy(self)
        other_organism.ncbi_tax_id = ncbi_tax_id
        return other_organism


class PickleMapping(object):
    
    def __init__(self, pickleFile):
        
        self.pickleFile = pickleFile


class ReadSettings:
    
    def __init__(
            self,
            name = "unknown",
            separator = None,
            nameColA = 0,
            nameColB = 1,
            name_typeA = "uniprot",
            name_typeB = "uniprot",
            typeA = "protein",
            typeB = "protein",
            isDirected = False,
            sign = False,
            inFile = None,
            references = False,
            extraEdgeAttrs = None,
            extraNodeAttrsA = None,
            extraNodeAttrsB = None,
            header = False,
            taxonA = 9606,
            taxonB = 9606,
            ncbiTaxId = False,
            interactionType = 'PPI',
            positiveFilters = None,
            negativeFilters = None,
            mark_source  =  None,
            mark_target  =  None,
            inputArgs = None,
            curlArgs = None,
            must_have_references = True,
            huge = False,
            resource = None,
            unique_fields = None
        ):
        """
        :param str mark_source:
            Creates a boolean vertex attribute and sets it True for the
            source vertex of directed interactions from this particular
            resource.
        :param str mark_target:
            Same as ``mark_source`` but for target vertices.
        """
        self.typeA = typeA
        self.typeB = typeB
        self.nameColA = nameColA
        self.nameColB = nameColB
        self.name_typeA = name_typeA
        self.name_typeB = name_typeB
        self.isDirected = isDirected
        self.inFile = inFile
        self.extraEdgeAttrs = extraEdgeAttrs or {}
        self.extraNodeAttrsA = extraNodeAttrsA or {}
        self.extraNodeAttrsB = extraNodeAttrsB or {}
        self.name = name
        self.separator = separator
        self.header = header
        self.refs = references
        self.sign = sign
        self.taxonA = taxonA
        self.taxonB = taxonB
        self.ncbiTaxId = ncbiTaxId
        self.intType = interactionType
        self.positiveFilters = positiveFilters or []
        self.negativeFilters = negativeFilters or []
        self.inputArgs = inputArgs or {}
        self.curlArgs = curlArgs or {}
        self.must_have_references = must_have_references and bool(references)
        self.huge = huge
        self.resource = self.name if resource is None else resource
        self.mark_source = mark_source
        self.mark_target = mark_target
        self.unique_fields = unique_fields or set()


class ReadList:
    def __init__(self,
                 name="unknown",
                 separator=None,
                 nameCol=0,
                 name_type="uniprot",
                 typ="protein",
                 inFile=None,
                 extraAttrs={},
                 header=False):
        self.typ = typ
        self.nameCol = nameCol
        self.name_type = name_type
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
