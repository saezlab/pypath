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
        self.input_args = input_args or {}
        self.ncbi_tax_id = ncbi_tax_id
        self.mapper = mapper
        
        self.load()
    
    
    def load(self):
        """
        Loads and preprocesses annotation data.
        """
        
        self.set_mapper()
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
        
        self._process_method()
    
    
    def load_uniprots(self):
        """
        Retrieves a set of all UniProt IDs to have a base set of the entire
        proteome.
        """
        
        self.uniprots = set(dataio.all_uniprots(organism = self.ncbi_tax_id))
    
    
    def _process_method(*args, **kwargs):
        
        pass
    
    
    def __contains__(self, uniprot):
        
        return uniprot in self.annot
    
    
    def __getitem__(self, uniprot):
        
        if uniprot in self:
            
            return self.annot[uniprot]


class Membranome(AnnotationBase):
    
    
    def __init__(self):
        
        AnnotationBase.__init__(
            self,
            name = 'Membranome',
            input_method = 'get_membranome',
        )
    
    
    def _process_method(self):
        
        record = collections.namedtuple(
            'MembranomeAnnotation',
            ['membrane', 'side'],
        )
        
        _annot = collections.defaultdict(set)
        
        for a in self.annot:
            
            _annot[a[0]].add(record(a[1], a[2]))
        
        self.annot = dict(_annot)


class Exocarta(AnnotationBase):
    
    
    def __init__(self, ncbi_tax_id = 9606, mapper = None, **kwargs):
        
        if 'organism' not in kwargs:
            
            kwargs['organism'] = ncbi_tax_id
        
        if 'database' not in kwargs:
            
            kwargs['database'] = 'exocarta'
        
        AnnotationBase.__init__(
            self,
            name = kwargs['database'].capitalize(),
            ncbi_tax_id = ncbi_tax_id,
            input_method = '_get_exocarta_vesiclepedia',
            mapper = mapper,
            input_args = kwargs,
        )
    
    
    def _process_method(self):
        
        record = collections.namedtuple(
            '%sAnnotation' % self.name,
            ['pmid', 'tissue', 'vesicle'],
        )
        
        _annot = collections.defaultdict(set)
        
        for a in self.annot:
            
            uniprots = self.mapper.map_name(a[1], 'genesymbol', 'uniprot')
            
            for u in uniprots:
                
                for vesicle in (
                    a[3][3] if self.name == 'Vesiclepedia' else ('Exosomes',)
                ):
                
                    _annot[u].add(record(a[3][0], a[3][2], vesicle))
        
        self.annot = dict(_annot)


class Vesiclepedia(Exocarta):
    
    
    def __init__(self, ncbi_tax_id = 9606, mapper = None, **kwargs):
        
        Exocarta.__init__(
            self,
            ncbi_tax_id = ncbi_tax_id,
            database = 'vesiclepedia',
            mapper = mapper,
            **kwargs,
        )


class Matrisome(AnnotationBase):
    
    
    def __init__(self, ncbi_tax_id = 9606, mapper = None, **kwargs):
        
        if 'organism' not in kwargs:
            
            kwargs['organism'] = ncbi_tax_id
        
        AnnotationBase.__init__(
            self,
            name = 'Matrisome',
            ncbi_tax_id = ncbi_tax_id,
            input_method = 'get_matrisome',
            input_args = kwargs,
            mapper = mapper,
        )
    
    
    def _process_method(self):
        
        _annot = collections.defaultdict(set)
        
        record = collections.namedtuple(
            'MatrisomeAnnotation',
            ['mainclass', 'subclass', 'subsubclass'],
        )
        
        for uniprot, a in iteritems(self.annot):
            
            _annot[uniprot].add(record(*a))
        
        self.annot = dict(_annot)


class Surfaceome(AnnotationBase):
    
    
    def __init__(self, mapper = None):
        
        AnnotationBase.__init__(
            self,
            name = 'Surfaceome',
            input_method = 'get_surfaceome',
            mapper = mapper,
        )
    
    
    def _process_method(self):
        
        _annot = collections.defaultdict(set)
        
        record = collections.namedtuple(
            'SurfaceomeAnnotation',
            ['score', 'mainclass', 'subclasses']
        )
        record.__defaults__ = (None, None)
        
        for uniprot, a in iteritems(self.annot):
            
            _annot[uniprot].add(
                record(
                    a[0],
                    a[1],
                    tuple(sorted(a[2])) if a[2] else None,
                )
            )
        
        self.annot = dict(_annot)


class CellSurfaceProteinAtlas(AnnotationBase):
    
    
    def __init__(self, ncbi_tax_id = 9606, mapper = None, **kwargs):
        """
        The name of this resource abbreviated as `CSPA`.
        """
        
        if 'organism' not in kwargs:
            
            kwargs['organism'] = ncbi_tax_id
        
        AnnotationBase.__init__(
            self,
            name = 'CSPA',
            ncbi_tax_id = ncbi_tax_id,
            input_method = 'get_cspa',
            process_method = self._process_method,
            input_args = kwargs,
            mapper = mapper,
        )
    
    
    @staticmethod
    def _process_method(annot, mapper):
        
        return dict((u, set()) for u in annot)


class HumanPlasmaMembraneReceptome(AnnotationBase):
    
    
    def __init__(self, mapper = None):
        """
        The name of this resource abbreviated as `HPMR`.
        """
        
        AnnotationBase.__init__(
            self,
            name = 'HPMR',
            input_method = 'get_hpmr',
            process_method = self._process_method,
            mapper = mapper,
        )
    
    
    def process_method(self):
        
        return dict(
            (uniprot, set())
            for genesymbol in annot
            for uniprot in mapper.map_name(
                genesymbol, 'genesymbol', 'uniprot'
            )
        )
