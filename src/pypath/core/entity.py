#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#  Exports tables for webservice.
#
#  Copyright (c) 2014-2020 - EMBL
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#                  Nicolàs Palacio
#                  Olga Ivanova
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

"""
Provides classes for representing molecular entities and their collections.
A molecular entity is defined by its identifier, type and taxon.
"""

from future.utils import iteritems

import importlib as imp
import collections

import pypath.share.common as common
import pypath.internals.intera as intera
import pypath.share.session as session_mod
import pypath.utils.mapping as mapping


EntityKey = collections.namedtuple(
    'EntityKey',
    [
        'identifier',
        'id_type',
        'entity_type',
        'taxon',
    ]
)


class Entity(session_mod.Logger):
    """
    Represents a molecular entity such as protein, miRNA, lncRNA or small
    molecule.
    
    :arg str identifier:
        An identifier from the reference database e.g. UniProt ID for
        proteins.
    :arg str entity_type:
        The type of the molecular entity, defaults to ``'protein'``.
    :arg str id_type:
        The type of the identifier (the reference database), default is
        ``'uniprot'``.
    :arg int taxon:
        The NCBI Taxonomy Identifier of the molecular entity, e.g. ``9606``
        for human. Use ``0`` for non taxon specific molecules e.g. metabolites
        or drug compounds.
    :arg NoneType,dict attrs:
        A dictionary of additional attributes.
    """
    
    __slots__ = [
        'identifier',
        'entity_type',
        'id_type',
        'taxon',
        'attrs',
        'label',
        'key',
    ]
    
    
    def __init__(
            self,
            identifier,
            entity_type = 'protein',
            id_type = 'uniprot',
            taxon = 9606,
            attrs = None,
        ):
        
        self.identifier = identifier
        self.entity_type = entity_type or self.get_entity_type()
        # override `protein` in case this is a `complex`
        self.entity_type = 'complex' if self.is_complex() else entity_type
        self.id_type = id_type
        self.taxon = taxon
        self.key = self._key
        
        self.attrs = attrs or {}
        
        self.set_label()
    
    
    def reload(self):
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        import importlib as imp
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    
    @staticmethod
    def entity_name_str(entity):
        
        return (
            entity
                if isinstance(entity, common.basestring) else
            str(entity)
        )
    
    
    @classmethod
    def igraph_vertex_name(cls, igraph_v):
        
        return cls.entity_name_str(igraph_v['name'])
    
    
    @staticmethod
    def igraph_vertex_label(igraph_v):
        
        return igraph_v['label']
    
    
    @classmethod
    def igraph_vertex_name_label(cls, igraph_v):
        
        return (
            cls.igraph_vertex_name(igraph_v),
            cls.igraph_vertex_label(igraph_v),
        )
    
    
    @staticmethod
    def _is_protein(key):
        
        return (
            isinstance(key, common.basestring) and
            not key.startswith('MIMAT') and
            not key.startswith('COMPLEX')
        )
    
    
    @staticmethod
    def _is_mirna(key):
        
        return (
            isinstance(key, common.basestring) and
            key.startswith('MIMAT')
        )
    
    
    @staticmethod
    def _is_complex(key):
        
        return isinstance(key, intera.Complex) or (
            isinstance(key, common.basestring) and
            key.startswith('COMPLEX')
        )
    
    
    @classmethod
    def _get_entity_type(cls, key):
        
        return (
            'complex'
                if cls._is_complex(key) else
            'mirna'
                if cls._is_mirna(key) else
            'protein'
        )
    
    
    def is_protein(self):
        
        return self._is_protein(self.identifier)
    
    
    def is_mirna(self):
        
        return self._is_mirna(self.identifier)
    
    
    def is_complex(self):
        
        return self._is_complex(self.identifier)
    
    
    def get_entity_type(self):
        
        return self._get_entity_type(self.identifier)
    
    
    @property
    def _key(self):
        
        return EntityKey(
            identifier = self.identifier,
            id_type = self.id_type,
            entity_type = self.entity_type,
            taxon = self.taxon,
        )
    
    
    def __hash__(self):
        
        return hash(self.key)
    
    
    def __eq__(self, other):
        
        return (
            self.__hash__() == other.__hash__()
                if hasattr(other, 'key') else
            self.identifier == other
        )
    
    
    def __lt__(self, other):
        
        return (
            self.key < other.key
                if hasattr(other, 'key') else
            self.identifier < other
        )
    
    
    def __gt__(self, other):
        
        return (
            self.key < other.key
                if hasattr(other, 'key') else
            self.identifier < other
        )
    
    
    def set_label(self):
        
        self.label = mapping.label(
            name = self.identifier,
            id_type = self.id_type,
            ncbi_tax_id = self.taxon,
        ) or self.identifier
    
    
    def __repr__(self):
        
        return '<Entity: %s>' % (self.label or self.identifier)
    
    
    def __iadd__(self, other):
        
        if self == other:
            
            self.update_attrs(**other.attrs)
        
        return self
    
    
    def update_attrs(self, **kwargs):
        
        for key, val in iteritems(kwargs):
            
            if key in self.attrs:
                
                self.attrs[key] = common.combine_attrs((self.attrs[key], val))
                
            else:
                
                self.attrs[key] = val


class EntityList(object):
    
    
    def __init__(self, entities):
        
        self._entities = (
            entities
                if isinstance(entities, (list, tuple, set)) else
            list(entities)
        )
    
    
    def __iter__(self):
        
        for e in self._entities:
            
            yield e
    
    
    def __len__(self):
        
        return len(self._entities)
    
    
    def __repr__(self):
        
        return '<EntityList (%u elements)>' % len(self)
    
    
    def __add__(self, other):
        
        return EntityList(set(itertools.chain(self._entities, list(other))))
    
    
    def __iadd__(self, other):
        
        self._entities = set(itertools.chain(self._entities, list(other)))
        
        return self
    
    
    @property
    def labels(self):
        
        for e in self:
            
            yield e.label
    
    
    @property
    def ids(self):
        
        for e in self:
            
            yield e.identifier
    
    
    identifiers = ids
    
    
    @property
    def entities(self):
        
        for e in self:
            
            yield e
    
    
    @property
    def list_ids(self):
        
        return list(self.ids)
    
    
    @property
    def list_labels(self):
        
        return list(self.labels)
    
    
    @property
    def list_entities(self):
        
        return list(self.entities)
    
    
    l = labels
    i = ids
    e = entities
    ll = list_labels
    li = list_ids
    le = list_entities
