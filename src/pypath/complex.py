#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2019
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#                  Nicolàs Palacio
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

import collections

import pypath.dataio as dataio
import pypath.intera as intera
import pypath.resource as resource


class AbstractComplexResource(resource.AbstractResource):
    """
    A resource which provides information about molecular complexes.
    """
    
    
    def __init__(
            self,
            name,
            ncbi_tax_id = 9606,
            input_method = None,
            input_args = None,
            **kwargs,
        ):
        """
        name : str
            Custom name for the resource.
        input_method : callable
            Method providing the input data.
        process_method : callable
            Method processing the data and yielding ``intera.Complex``
            instances.
        """
        
        self.complexes = {}
        
        resource.AbstractResource.__init__(
            self,
            name = name,
            ncbi_tax_id = ncbi_tax_id,
            input_method = input_method,
            input_args = input_args,
        )
        
        self.load()
    
    
    def load(self):
        
        resource.AbstractResource.load(self)
        self.update_index()
    
    
    def _process_method(self):
        
        self.complexes = self.data
        
        delattr(self, 'data')
    
    
    def __iter__(self):
        
        for cplex in self.complexes.values():
            
            yield cplex
    
    
    def update_index(self):
        
        self.proteins = collections.defaultdict(set)
        self.resources = collections.defaultdict(set)
        
        for cplex in self:
            
            for protein in cplex:
                
                self.proteins[protein].add(cplex)
            
            for db in cplex.sources:
                
                self.resources[protein].add(cplex)
    
    
    def __contains__(self, other):
        
        # a Complex instance
        if isinstance(other, intera.Complex):
            
            other = other.__str__()
        
        # either a UniProt ID or
        # a complex string representation
        if isinstance(other, common.basestring):
            
            if len(other) <= 10:
                
                return other in self.proteins
                
            else:
                
                return other in self.complexes
        
        return False


class CellPhoneDB(AbstractComplexResource):
    
    
    def __init__(self, **kwargs):
        
        AbstractComplexResource.__init__(
            self,
            name = 'CellPhoneDB',
            input_method = 'cellphonedb_complexes',
        )


class Corum(AbstractComplexResource):
    
    
    def __init__(self, input_args = None, **kwargs):
        
        AbstractComplexResource.__init__(
            self,
            name = 'CORUM',
            input_method = 'corum_complexes',
            input_args = input_args or {},
        )


class Havugimana(AbstractComplexResource):
    
    
    def __init__(self, input_args = None, **kwargs):
        
        AbstractComplexResource.__init__(
            self,
            name = 'Havugimana2012',
            input_method = 'havugimana_complexes',
            input_args = input_args or {},
        )


class Compleat(AbstractComplexResource):
    
    
    def __init__(self, input_args = None, **kwargs):
        
        AbstractComplexResource.__init__(
            self,
            name = 'Compleat',
            input_method = 'compleat_complexes',
            input_args = input_args or {},
        )


class ComplexPortal(AbstractComplexResource):
    
    
    def __init__(self, input_args = None, **kwargs):
        
        AbstractComplexResource.__init__(
            self,
            name = 'ComplexPortal',
            input_method = 'complexportal_complexes',
            input_args = input_args or {},
        )


class Pdb(AbstractComplexResource):
    
    
    def __init__(self, input_args = None, **kwargs):
        
        AbstractComplexResource.__init__(
            self,
            name = 'PDB',
            input_method = 'pdb_complexes',
            input_args = input_args or {},
        )


class Signor(AbstractComplexResource):
    
    
    def __init__(self, input_args = None, **kwargs):
        
        AbstractComplexResource.__init__(
            self,
            name = 'Signor',
            input_method = 'signor_complexes',
            input_args = input_args or {},
        )
