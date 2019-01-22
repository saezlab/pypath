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
from past.builtins import xrange, range, reduce


import imp
import collections


import pypath.dataio as dataio
import pypath.common as common
import pypath.mapping as mapping
import pypath.go as go
import pypath.intercell_annot as intercell_annot


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
        Calls the ``_process_method``.
        """
        
        self._process_method()
    
    
    def load_uniprots(self):
        """
        Retrieves a set of all UniProt IDs to have a base set of the entire
        proteome.
        """
        
        self.uniprots = set(dataio.all_uniprots(organism = self.ncbi_tax_id))
    
    
    def _process_method(self, *args, **kwargs):
        """
        By default it converts a set to dict of empty sets in order to make
        it compatible with other methods.
        Derived classes might override.
        """
        
        self.annot = dict((u, set()) for u in self.annot)
    
    
    def get_subset(self, **kwargs):
        """
        Retrieves a subset by filtering based on ``kwargs``.
        Each argument should be a name and a value or set of values.
        Elements having the provided values in the annotation will be
        returned.
        Returns a set of UniProt IDs.
        """
        
        result = set()
        
        for uniprot, annot in iteritems(self.annot):
            
            for name, value in iteritems(kwargs):
                
                if not any(
                    (
                        getattr(a, name) == value or (
                            isinstance(
                                getattr(a, name),
                                (common.basestring, tuple)
                            ) and
                            getattr(a, name) in value
                        )
                    )
                    for a in annot
                ):
                    
                    break
                
                result.add(uniprot)
        
        return result
    
    
    def to_set(self):
        
        return set(self.annot.keys())
    
    
    def coverage(self, other):
        
        other = other if isinstance(other, set) else set(other)
        
        return len(self & other) / len(self)
    
    
    def proportion(self, other):
        
        other = other if isinstance(other, set) else set(other)
        
        return len(self & other) / len(other)
    
    
    def subset_intersection(self, universe, **kwargs):
        
        subset = self.get_subset(**kwargs)
        
        return len(subset & universe) / len(subset)
    
    
    def get_values(self, name, exclude_none = True):
        
        values =  set(
            getattr(a, name)
            for aset in self.annot.values()
            for a in aset
        )
        
        if exclude_none:
            
            values.discard(None)
        
        return values
    
    
    def get_names(self):
        
        names = ()
        
        for val in self.annot.values():
            
            if val:
                
                names = val._fields
                break
        
        return names
    
    
    def __contains__(self, uniprot):
        
        return uniprot in self.annot
    
    
    def __getitem__(self, uniprot):
        
        if uniprot in self:
            
            return self.annot[uniprot]
    
    
    def __and__(self, other):
        
        return other & self.to_set()
    
    
    def __or__(self, other):
        
        return other | self.to_set()
    
    
    def __sub__(self, other):
        
        return self.to_set() - other
    
    
    def __len__(self):
        
        return len(self.annot)


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
            input_args = kwargs,
            mapper = mapper,
        )


class HumanPlasmaMembraneReceptome(AnnotationBase):
    
    
    def __init__(self, mapper = None):
        """
        The name of this resource abbreviated as `HPMR`.
        """
        
        AnnotationBase.__init__(
            self,
            name = 'HPMR',
            input_method = 'get_hpmr',
            mapper = mapper,
        )
    
    
    def _process_method(self):
        
        self.annot = dict(
            (uniprot, set())
            for genesymbol in self.annot
            for uniprot in self.mapper.map_name(
                genesymbol, 'genesymbol', 'uniprot'
            )
        )


class MatrixdbBase(AnnotationBase):
    
    
    def __init__(self, category, ncbi_tax_id = 9606):
        """
        Protein annotations from MatrixDB.
        
        :arg str category:
            The protein annotation category. Possible values: `ecm`, `membrane`
            or `secreted`.
        """
        
        AnnotationBase.__init__(
            self,
            name = 'MatrixDB_%s' % category,
            ncbi_tax_id = ncbi_tax_id,
            input_method = 'matrixdb_%s_proteins' % category.lower(),
        )


class MatrixdbSecreted(MatrixdbBase):
    
    
    def __init__(self, ncbi_tax_id = 9606):
        """
        Secreted proteins annotations from MatrixDB.
        
        :arg int ncbi_tax_id:
            NCBI Taxonomy ID of the organism.
        """
        
        MatrixdbAnnotation.__init__(
            self,
            category = 'Secreted',
            ncbi_tax_id = ncbi_tax_id,
        )


class MatrixdbMembrane(MatrixdbBase):
    
    
    def __init__(self, ncbi_tax_id = 9606):
        """
        Membrane proteins annotations from MatrixDB.
        
        :arg int ncbi_tax_id:
            NCBI Taxonomy ID of the organism.
        """
        
        MatrixdbAnnotation.__init__(
            self,
            category = 'Membrane',
            ncbi_tax_id = ncbi_tax_id,
        )


class MatrixdbECM(MatrixdbBase):
    
    
    def __init__(self, ncbi_tax_id = 9606):
        """
        ECM proteins annotations from MatrixDB.
        
        :arg int ncbi_tax_id:
            NCBI Taxonomy ID of the organism.
        """
        
        MatrixdbAnnotation.__init__(
            self,
            category = 'ECM',
            ncbi_tax_id = ncbi_tax_id,
        )


class Locate(AnnotationBase):
    
    
    def __init__(
            self,
            ncbi_tax_id = 9606,
            mapper = None,
            literature = True,
            external = True,
            predictions = False,
        ):
        
        input_args = {
            'organism': ncbi_tax_id or 9606,
            'mapper': mapper,
            'literature': literature,
            'external': external,
            'predictions': predictions,
        }
        
        AnnotationBase.__init__(
            self,
            name = 'Locate',
            input_method = 'get_locate_localizations',
            ncbi_tax_id = ncbi_tax_id,
            input_args = input_args,
        )
    
    
    def _process_method(self, *args, **kwargs):
        
        #  already the appropriate format, no processing needed
        pass


class GOCustomIntercell(go.GOCustomAnnotation):
    
    
    def __init__(
            self,
            categories = None,
            go_annot = None,
            ncbi_tax_id = 9606,
        ):
        """
        Same as :class:``pypath.go.GOCustomAnnotation``
        initialized with the categories defined in
        ``pypath.intercell_annot.intercell_categories``.
        """
        
        categories = categories or intercell_annot.intercell_categories
        
        go.GOCustomAnnotation.__init__(
            self,
            categories = categories,
            go_annot = go_annot,
            ncbi_tax_id = ncbi_tax_id,
        )


class GOIntercell(AnnotationBase):
    
    
    def __init__(
            self,
            categories = None,
            go_annot = None,
            ncbi_tax_id = 9606,
        ):
        """
        Annotation of proteins based on their roles in intercellular
        communication from Gene Ontology.
        """
        
        self.categories = categories
        self.go_annot = go_annot
        
        AnnotationBase.__init__(
            self,
            name = 'GO_Intercell',
            ncbi_tax_id = ncbi_tax_id,
        )
    
    
    def load(self):
        
        record = collections.namedtuple(
            'GOIntercellAnnotation',
            ('mainclass',),
        )
        
        annot = GOCustomIntercell(
            categories = self.categories,
            go_annot = self.go_annot,
            ncbi_tax_id = self.ncbi_tax_id,
        )
        
        annot_uniprots = annot.get_annotations()
        
        _annot = collections.defaultdict(set)
        
        for mainclass, uniprots in iteritems(annot_uniprots):
            
            for uniprot in uniprots:
                
                _annot[uniprot].add(record(mainclass = mainclass))
        
        self.annot = dict(_annot)
    
    
    def _process_method(self, *args, **kwargs):
        
        pass
