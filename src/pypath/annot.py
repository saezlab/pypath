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


import sys
import imp
import collections
import itertools

import numpy as np
import pandas as pd

import pypath.dataio as dataio
import pypath.common as common
import pypath.mapping as mapping
import pypath.resource as resource
import pypath.go as go
import pypath.intercell_annot as intercell_annot


annotation_sources = {
    'Membranome',
    'Exocarta',
    'Vesiclepedia',
    'Matrisome',
    'Surfaceome',
    'CellSurfaceProteinAtlas',
    'HumanPlasmaMembraneReceptome',
    'MatrixdbSecreted',
    'MatrixdbMembrane',
    'MatrixdbECM',
    'Locate',
    'GOIntercell',
    'CellPhoneDB',
}

default_fields = {
    'Matrisome': ('mainclass', 'subclass'),
    'Locate': ('location',),
    'Vesiclepedia': ('vesicle',),
    'Exocarta': ('vesicle',),
    'CellPhoneDB': (
        'receptor',
        'adhesion',
        'cytoplasm',
        'peripheral',
        'secretion',
        'secreted',
        'transporter',
        'transmembrane',
        'extracellular',
    )
}


class AnnotationBase(resource.AbstractResource):
    
    
    def __init__(
            self,
            name,
            ncbi_tax_id = 9606,
            input_method = None,
            input_args = None,
            **kwargs,
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
        
        resource.AbstractResource.__init__(
            self,
            name = name,
            ncbi_tax_id = ncbi_tax_id,
            input_method = input_method,
            input_args = input_args,
        )
        
        self.load()
    
    
    def reload(self):
        """
        Reloads the object from the module level.
        """
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    
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
        
        self.annot = dict((u, set()) for u in self.data)
    
    
    def get_subset(self, method = None, **kwargs):
        """
        Retrieves a subset by filtering based on ``kwargs``.
        Each argument should be a name and a value or set of values.
        Elements having the provided values in the annotation will be
        returned.
        Returns a set of UniProt IDs.
        """
        
        result = set()
        
        for uniprot, annot in iteritems(self.annot):
            
            for a in annot:
                
                # we either call a method on all records
                # or check against conditions provided in **kwargs
                if (
                    callable(method) and
                    method(a)
                ) or all(
                    (
                        getattr(a, name) == value or (
                            callable(value) and
                            value(getattr(a, name))
                        ) or (
                            isinstance(
                                getattr(a, name),
                                (common.basestring, tuple)
                            ) and
                            getattr(a, name) in value
                        )
                    )
                    for name, value in iteritems(kwargs)
                ):
                    
                    result.add(uniprot)
                    break
        
        return result
    
    
    def get_subset_bool_array(self, uniprots, **kwargs):
        
        subset = self.get_subset(**kwargs)
        
        return np.array([
            uniprot in subset
            for uniprot in uniprots
        ])
    
    
    def to_bool_array(self, uniprots):
        
        total = self.to_set()
        
        return np.array([
            uniprot in total
            for uniprot in uniprots
        ])
    
    
    def to_set(self):
        
        return set(self.annot.keys())
    
    
    def all_uniprots(self):
        """
        All UniProt IDs annotated in this resource.
        """
        
        return sorted(self.annot.keys())
    
    
    def to_array(self, uniprots = None, use_fields = None):
        
        uniprots = uniprots or self.all_uniprots()
        all_fields = self.get_names()
        fields = use_fields or all_fields
        ifields = tuple(
            i for i, field in enumerate(all_fields) if field in fields
        )
        result = [
            (
                (self.name,),
                self.to_bool_array(uniprots = uniprots)
            )
        ]
        
        for i in xrange(len(fields)):
            
            this_ifields = ifields[:i+1]
            this_fields  =  fields[:i+1]
            
            value_combinations = set(
                tuple(annot[j] for j in this_ifields)
                for annots in self.annot.values()
                for annot in annots
            )
            value_combinations = sorted(
                values
                for values in value_combinations
                if not any(v is None for v in values) and
                not any(isinstance(v, float) for v in values)
            )
            
            for values in value_combinations:
                
                this_values = dict(zip(this_fields, values))
                
                this_array = self.get_subset_bool_array(
                    uniprots = uniprots,
                    **this_values,
                )
                
                result.append(
                    (
                        (self.name,) + values,
                        this_array,
                    )
                )
        
        return (
            tuple(r[0] for r in result),
            np.vstack([r[1] for r in result]).T
        )
    
    
    def make_df(self):
        
        discard = {'n/a', None}
        
        columns = [
            'uniprot',
            'source',
            'label',
            'value',
            'record_id',
        ]
        
        records = []
        
        irec = 0
        
        for uniprot, annots in iteritems(self.annot):
            
            if not annots:
                
                records.append([
                    uniprot,
                    self.name,
                    'in %s' % self.name,
                    'yes',
                    irec,
                ])
                
                irec += 1
            
            for annot in annots:
                
                for label, value in zip(annot._fields, annot):
                    
                    if value in discard:
                        
                        continue
                    
                    if isinstance(value, (set, list, tuple)):
                        
                        value = ';'.join(map(str, value))
                    
                    records.append([
                        uniprot,
                        self.name,
                        label,
                        str(value),
                        irec,
                    ])
                
                irec += 1
        
        self.df = pd.DataFrame(
            records,
            columns = columns,
        )
    
    
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
        
        for values in self.annot.values():
            
            if values:
                
                for val in values:
                    
                    names = val._fields
                    break
            
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
    
    
    def __init__(self, **kwargs):
        
        AnnotationBase.__init__(
            self,
            name = 'Membranome',
            input_method = 'get_membranome',
            **kwargs,
        )
    
    
    def _process_method(self):
        
        record = collections.namedtuple(
            'MembranomeAnnotation',
            ['membrane', 'side'],
        )
        
        _annot = collections.defaultdict(set)
        
        for a in self.data:
            
            _annot[a[0]].add(record(a[1], a[2]))
        
        self.annot = dict(_annot)


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
            input_args = kwargs,
        )
    
    
    def _process_method(self):
        
        record = collections.namedtuple(
            '%sAnnotation' % self.name,
            ['pmid', 'tissue', 'vesicle'],
        )
        
        _annot = collections.defaultdict(set)
        
        for a in self.data:
            
            uniprots = mapping.map_name(a[1], 'genesymbol', 'uniprot')
            
            for u in uniprots:
                
                for vesicle in (
                    a[3][3] if self.name == 'Vesiclepedia' else ('Exosomes',)
                ):
                
                    _annot[u].add(record(a[3][0], a[3][2], vesicle))
        
        self.annot = dict(_annot)


class Vesiclepedia(Exocarta):
    
    
    def __init__(self, ncbi_tax_id = 9606, **kwargs):
        
        Exocarta.__init__(
            self,
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
    
    
    def _process_method(self):
        
        _annot = collections.defaultdict(set)
        
        record = collections.namedtuple(
            'MatrisomeAnnotation',
            ['mainclass', 'subclass', 'subsubclass'],
        )
        
        for uniprot, a in iteritems(self.data):
            
            _annot[uniprot].add(record(*a))
        
        self.annot = dict(_annot)


class Surfaceome(AnnotationBase):
    
    
    def __init__(self, **kwargs):
        
        AnnotationBase.__init__(
            self,
            name = 'Surfaceome',
            input_method = 'get_surfaceome',
            **kwargs,
        )
    
    
    def _process_method(self):
        
        _annot = collections.defaultdict(set)
        
        record = collections.namedtuple(
            'SurfaceomeAnnotation',
            ['score', 'mainclass', 'subclasses']
        )
        record.__defaults__ = (None, None)
        
        for uniprot, a in iteritems(self.data):
            
            _annot[uniprot].add(
                record(
                    a[0],
                    a[1],
                    tuple(sorted(a[2])) if a[2] else None,
                )
            )
        
        self.annot = dict(_annot)


class CellSurfaceProteinAtlas(AnnotationBase):
    
    
    def __init__(self, ncbi_tax_id = 9606, **kwargs):
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
        )


class HumanPlasmaMembraneReceptome(AnnotationBase):
    
    
    def __init__(self, **kwargs):
        """
        The name of this resource abbreviated as `HPMR`.
        """
        
        AnnotationBase.__init__(
            self,
            name = 'HPMR',
            input_method = 'get_hpmr',
            **kwargs,
        )
    
    
    def _process_method(self):
        
        self.annot = dict(
            (uniprot, set())
            for genesymbol in self.data
            for uniprot in mapping.map_name(
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
        
        MatrixdbBase.__init__(
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
        
        MatrixdbBase.__init__(
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
        
        MatrixdbBase.__init__(
            self,
            category = 'ECM',
            ncbi_tax_id = ncbi_tax_id,
        )


class Locate(AnnotationBase):
    
    
    def __init__(
            self,
            ncbi_tax_id = 9606,
            literature = True,
            external = True,
            predictions = False,
        ):
        
        input_args = {
            'organism': ncbi_tax_id or 9606,
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
    
    
    def _process_method(self):
        
        #  already the appropriate format, no processing needed
        self.annot = self.data
        
        delattr(self, 'data')


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


class CellPhoneDB(AnnotationBase):
    
    
    def __init__(self, **kwargs):
        
        AnnotationBase.__init__(
            self,
            name = 'CellPhoneDB',
            input_method = 'cellphonedb_protein_annotations',
            ncbi_tax_id = 9606,
        )
    
    
    def _process_method(self, *args, **kwargs):
        
        self.annot = dict(
            (uniprot, {annot,})
            for uniprot, annot in
            iteritems(self.data)
        )


class AnnotationTable(object):
    
    
    def __init__(
            self,
            uniprots = None,
            use_sources = None,
            use_fields = None,
            ncbi_tax_id = 9606,
            swissprot_only = True,
            keep_annotators = False,
            load = True,
        ):
        """
        Sorry Nico I don't write docs because lab meeting tomorrow!
        """
        
        self._module = sys.modules[self.__module__]
        self.use_sources = use_sources or annotation_sources
        self.use_fields = use_fields or default_fields
        self.ncbi_tax_id = ncbi_tax_id
        self.keep_annotators = keep_annotators
        self.uniprots = (
            uniprots or
            sorted(
                dataio.all_uniprots(
                    organism = ncbi_tax_id,
                    swissprot = swissprot_only,
                )
            )
        )
        self.rows = dict(
            reversed(i)
            for i in enumerate(self.uniprots)
        )
        
        if load:
            
            self.load()
    
    
    def reload(self):
        """
        Reloads the object from the module level.
        """
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    
    def load(self):
        
        annots = {}
        names  = []
        arrays = []
        
        for cls in self.use_sources:
            
            annot = getattr(self._module, cls)(
                ncbi_tax_id = self.ncbi_tax_id
            )
            
            use_fields = (
                self.use_fields[cls] if cls in self.use_fields else None
            )
            
            this_names, this_array = annot.to_array(
                uniprots = self.uniprots,
                use_fields = use_fields
            )
            
            names.extend(this_names)
            arrays.append(this_array)
            
            if self.keep_annotators:
                
                annots[annot.name] = annot
        
        self.annots = annots
        self.names = np.array(list(itertools.chain(names)))
        self.data = np.hstack(arrays)
        self.set_cols()
        self.uniprots = np.array(self.uniprots)
    
    
    def set_cols(self):
        
        self.cols = dict((name, i) for i, name in enumerate(self.names))
    
    
    def keep(self, keep):
        
        ikeep = np.array([
            i for i, name in enumerate(self.names) if name in keep
        ])
        
        self.names = self.names[ikeep]
        self.data  = self.data[:,ikeep]
        self.set_cols()
    
    
    def make_sets(self):
        
        self.sets = dict(
            (
                name,
                set(self.uniprots[self.data[:,i]])
            )
            for i, name in enumerate(self.names)
        )
    
    
    def annotate_network(self, pa):
        
        nodes = pa.graph.vs['name']
        edges = [
            (
                nodes[e.source],
                nodes[e.target]
            )
            for e in pa.graph.es
        ]
        
        nodeannot = []
        edgeannot = []
        
        for i, uniprot in enumerate(nodes):
            
            for name, uniprots in iteritems(self.sets):
                
                if uniprot in uniprots:
                    
                    nodeannot.append((name, i))
        
        for i, (uniprot1, uniprot2) in enumerate(edges):
            
            for name1, uniprots1 in iteritems(self.sets):
                
                for name2, uniprots2 in iteritems(self.sets):
                    
                    if uniprot1 in uniprots1 and uniprot2 in uniprots2:
                        
                        edgeannot.append((name1, name2, i))
        
        return nodeannot, edgeannot
    
    
    def network_stats(self, pa):
        
        nodeannot, edgeannot = self.annotate_network(pa)
        
        nodestats = collections.Counter('__'.join(n[0]) for n in nodeannot)
        
        edgestats = collections.Counter(
            tuple(sorted(('__'.join(e[0]), '__'.join(e[1]))))
            for e in edgeannot
        )
        
        return nodestats, edgestats
    
    
    def export_network_stats(self, pa):
        
        nodestats, edgestats = self.network_stats(pa)
        
        with open('annot_edgestats2.tsv', 'w') as fp:
            
            _ = fp.write('\t'.join(('name1', 'name2', 'count')))
            _ = fp.write('\n')
            
            _ = fp.write('\n'.join(
                '%s\t%s\t%u' % (name1, name2, cnt)
                for (name1, name2), cnt in iteritems(edgestats)
            ))
        
        with open('annot_nodestats2.tsv', 'w') as fp:
            
            _ = fp.write('\t'.join(('name', 'count')))
            _ = fp.write('\n')
            
            _ = fp.write('\n'.join(
                '%s\t%u' % (name, cnt)
                for name, cnt in iteritems(nodestats)
            ))
    
    
    def to_dataframe(self):
        
        colnames = ['__'.join(name) for name in self.names]
        
        df = pd.DataFrame(
            data = self.data,
            index = self.uniprots,
            columns = colnames,
        )
        
        return df
    
    
    def make_narrow_df(self):
        
        for annot in self.annots.values():
            
            annot.make_df()
        
        self.narrow_df = pd.concat(
            annot.df for annot in self.annots.values()
        )


def init_db(keep_annotators = True):
    """
    Initializes or reloads the annotation database.
    The database will be assigned to the ``db`` attribute of this module.
    """
    
    globals()['db'] = AnnotationTable(keep_annotators = keep_annotators)


def get_db(keep_annotators = True):
    """
    Retrieves the current database instance and initializes it if does
    not exist yet.
    """
    
    if 'db' not in globals():
        
        init_db(keep_annotators = keep_annotators)
    
    return globals()['db']
