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

from future.utils import iteritems

import collections

import numpy as np
import pandas as pd

import pypath.dataio as dataio
import pypath.intera as intera
import pypath.resource as resource
import pypath.settings as settings


complex_resources = (
    'Signor',
    'Corum',
    'CellPhoneDB',
    'Havugimana',
    'Compleat',
    'ComplexPortal',
    'Pdb',
    'Hpmr',
    'GuideToPharmacology',
)


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
    
    
    def make_df(self):
        
        have_stoichiometry = {
            'PDB',
            'Compleat',
            'ComplexPortal',
            'CellPhoneDB',
        }
        
        colnames = [
            'name',
            'all_components',
            'component_uniprot',
            'component_stoichiometry',
            'sources',
            'references',
        ]
        
        records = []
        
        for cplex in self.complexes.values():
            
            has_stoi = have_stoichiometry & cplex.sources
            
            for comp, stoi in iteritems(cplex.components):
                
                records.append([
                    cplex.name if cplex.name else None,
                    cplex.__str__(),
                    comp,
                    '%u' % stoi if has_stoi else np.nan,
                    ';'.join(cplex.sources),
                    ';'.join(cplex.references),
                ])
        
        self.df = pd.DataFrame(
            records,
            columns = colnames,
        )


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
        
        input_args = input_args or {}
        
        if 'organism' not in input_args:
            
            input_args['organism'] = settings.get('default_organism')
        
        AbstractComplexResource.__init__(
            self,
            name = 'PDB',
            input_method = 'pdb_complexes',
            input_args = input_args or {},
        )


class Signor(AbstractComplexResource):
    
    
    def __init__(self, input_args = None, **kwargs):
        
        input_args = input_args or {}
        
        if 'organism' not in input_args:
            
            input_args['organism'] = settings.get('default_organism')
        
        AbstractComplexResource.__init__(
            self,
            name = 'Signor',
            input_method = 'signor_complexes',
            input_args = input_args or {},
        )


class Hpmr(AbstractComplexResource):
    
    
    def __init__(self, input_args = None, **kwargs):
        
        input_args = input_args or {}
        
        AbstractComplexResource.__init__(
            self,
            name = 'HPMR',
            input_method = 'hpmr_complexes',
            input_args = input_args or {},
        )


class GuideToPharmacology(AbstractComplexResource):
    
    
    def __init__(self, input_args = None, **kwargs):
        
        input_args = input_args or {}
        
        AbstractComplexResource.__init__(
            self,
            name = 'Guide2Pharma',
            input_method = 'guide2pharma_complexes',
            input_args = input_args or {},
        )


class ComplexAggregator(AbstractComplexResource):
    
    
    def __init__(
            self,
            resources = None,
        ):
        """
        Combines complexes from multiple resources.
        
        :arg list resources:
            List of resources. Names of complex resource classes in this
            module or custom 
        """
        
        self.resources = resources or complex_resources
        
        AbstractComplexResource.__init__(
            self,
            name = 'OmniPath',
        )
    
    
    def load(self):
        
        self.data = {}
        
        for res in self.resources:
            
            if not callable(res):
                
                if res in globals():
                    
                    res = globals()[res]
            
            if callable(res):
                
                processor = res()
                
            elif hasattr(res, 'complexes'):
                
                processor = res
            
            for key, cplex in iteritems(processor.complexes):
                
                if key in self.data:
                    
                    self.data[key] += cplex
                    
                else:
                    
                    self.data[key] = cplex
        
        resource.AbstractResource.load(self)
        self.update_index()


def init_db():
    """
    Initializes or reloads the complex database.
    The database will be assigned to the ``db`` attribute of this module.
    """
    
    globals()['db'] = ComplexAggregator()


def get_db():
    """
    Retrieves the current database instance and initializes it if does
    not exist yet.
    """
    
    if 'db' not in globals():
        
        init_db()
    
    return globals()['db']
