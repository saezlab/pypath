#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2018
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
from past.builtins import xrange, range, reduce


import imp
import collections


import pypath.dataio as dataio
import pypath.common as common
import pypath.mapping as mapping


class AnnotationBase(object):
    
    
    def __init__(
            self,
            name,
            mapper = None,
            ncbi_tax_id = 9606,
            input_method = None,
            process_method = None,
            input_args = None,
        ):
        """
        Represents annotations for a set of proteins.
        Loads the data from the original resource and provides methods
        to query the annotations.
        
        :arg str name:
            A custom name for the annotation resource.
        :arg int ncbi_tax_id:
            NCBI Taxonomy identifier.
        :arg callable,str input_method:
            Either a callable or the name of a method in the ``dataio``
            module. Should return a dict with UniProt IDs as keys or an
            object suitable for ``process_method``.
        :arg dict input_args:
            Arguments for the ``input_method``.
        """
        
        self.name = name
        self._input_method = input_method
        self._process_method = process_method or lambda x: x
        self.input_args = input_args or {}
        self.ncbi_tax_id = ncbi_tax_id
        self.mapper = mapper
        
        self.load()
    
    
    def load(self):
        """
        Loads and preprocesses annotation data.
        """
        
        self.set_mapper()
        self.load_uniprots()
        self.load_data()
        self.process()
    
    
    def set_mapper(self):
        
        if self.mapper is None:
            
            self.mapper = mapping.Mapper(ncbi_tax_id = self.ncbi_tax_id)
    
    
    def set_method(self):
        """
        Sets the data input method by looking up in ``dataio`` module if
        necessary.
        """
        
        if (
            isinstance(self._input_method, common.basestring) and
            hasattr(dataio, self._input_method)
        ):
            
            self.input_method = getattr(dataio, self._input_method)
    
    
    def load_data(self):
        """
        Loads the data by calling ``input_method``.
        """
        
        self.set_method()
        
        if hasattr(self, 'input_method'):
            
            self.annot = self.input_method(**self.input_args)
    
    
    def process(self):
        """
        Does nothing, derived classes might override.
        """
        
        self.annot = self.process_method(self.annot)
    
    
    def load_uniprots(self):
        """
        Retrieves a set of all UniProt IDs to have a base set of the entire
        proteome.
        """
        
        self.uniprots = set(datio.all_uniprots(organism = self.ncbi_tax_id))
    
    
    def __contains__(self, uniprot):
        
        return uniprot in self.annot
    
    
    def __getitem__(self, uniprot):
        
        if uniprot in self:
            
            return self.annot[uniprot]


class Membranome(AnnotationBase):
    
    
    def __init__(self, ncbi_tax_id = 9606, **kwargs):
        
        if 'organism' not in kwargs:
            
            kwargs['organism'] = ncbi_tax_id
        
        AnnotationBase.__init__(
            self,
            name = 'Membranome',
            ncbi_tax_id = ncbi_tax_id,
            input_method = 'get_membranome',
            process_method = self._process_method,
            input_args = kwargs,
        )
    
    
    @staticmethod
    def _process_method(annot):
        
        _annot = collections.defaultdict(set)
        
        for a in annot:
            
            _annot[a[0]].add(a[1:])
        
        return dict(_annot)


class Exocarta(AnnotationBase):
    
    
    def __init__(self, ncbi_tax_id = 9606, **kwargs):
        
        if 'organism' not in kwargs:
            
            kwargs['organism'] = ncbi_tax_id
        
        if 'database' not in kwargs:
            
            kwargs['database'] = 'exocarta'
        
        AnnotationBase.__init__(
            self,
            name = kwargs['database'].capitalize(),
            ncbi_tax_id = ncbi_tax_id,
            input_method = '_get_exocarta_vesiclepedia',
            process_method = self._process_method,
            input_args = kwargs,
        )
    
    
    @staticmethod
    def _process_method(annot):
        
        _annot = collections.defaultdict(set)
        
        for a in annot:
            
            uniprots = self.mapper.map_name(a[1], 'genesymbol', 'uniprot')
            
            for u in uniprots:
                
                _annot[u].add(a[3])
        
        return dict(_annot)


class Vesiclepedia(Exocarta):
    
    
    def __init__(self, ncbi_tax_id = 9606, **kwargs):
        
        Exocarta.__init__(
            ncbi_tax_id = ncbi_tax_id,
            database = 'vesiclepedia',
            **kwargs,
        )


class Matrisome(AnnotationBase):
    
    
    def __init__(self, ncbi_tax_id = 9606, **kwargs):
        
        if 'organism' not in kwargs:
            
            kwargs['organism'] = ncbi_tax_id
        
        AnnotationBase.__init__(
            self,
            name = 'Matrisome',
            ncbi_tax_id = ncbi_tax_id,
            input_method = 'get_matrisome',
            input_args = kwargs,
        )


class Surfaceome(AnntationBase):
    
    
    def __init__(self):
        
        if 'organism' not in kwargs:
            
            kwargs['organism'] = ncbi_tax_id
        
        AnnotationBase.__init__(
            self,
            name = 'Surfaceome',
            ncbi_tax_id = ncbi_tax_id,
            input_method = 'get_surfaceome',
        )
