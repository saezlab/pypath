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

import itertools
import importlib as imp
import collections

import pypath.share.common as common
from pypath.internals import intera
import pypath.share.session as session_mod
import pypath.utils.mapping as mapping
import pypath.utils.uniprot as uniprot
import pypath.share.settings as settings


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


    _default_id_types = settings.get('default_name_types')

    _id_type_to_entity_type = {
        'uniprot': 'protein',
        'genesymbol': 'protein',
        'mir-name': 'mirna',
        'mir-mat-name': 'mirna',
        'mir-pre': 'mirna',
        'mir-mat': 'mirna',
        'lncrna-genesymbol': 'lncrna',
    }

    _label_types = set(mapping.Mapper.label_type_to_id_type.keys())


    def __init__(
            self,
            identifier,
            entity_type = None,
            id_type = None,
            taxon = 9606,
            attrs = None,
        ):

        if (
            isinstance(identifier, Entity) or
            hasattr(identifier, 'identifier')
        ):

            (
                identifier,
                entity_type,
                id_type,
                taxon,
            ) = (
                identifier.identifier,
                identifier.entity_type,
                identifier.id_type,
                identifier.taxon,
            )

        self._bootstrap(identifier, id_type, entity_type, taxon)
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


    def _bootstrap(self, identifier, id_type, entity_type, taxon):

        if self._is_complex(identifier):

            entity_type = 'complex'
            id_type = 'complex'
            taxon = (
                identifier.ncbi_tax_id
                    if hasattr(identifier, 'ncbi_tax_id') else
                taxon
            )

        taxon = taxon or settings.get('default_organism')

        if not entity_type:

            if id_type and id_type in self._id_type_to_entity_type:

                entity_type = self._id_type_to_entity_type[id_type]


        if not id_type:

            id_type, entity_type = mapping.guess_type(
                identifier,
                entity_type = entity_type,
            )

        if not id_type and (not entity_type or entity_type == 'protein'):

            id_type, entity_type = 'genesymbol', 'protein'

        if id_type in self._label_types:

            _identifier = mapping.id_from_label0(
                label = identifier,
                label_id_type = id_type,
                ncbi_tax_id = taxon,
            )



            if _identifier and _identifier != identifier:

                id_type = mapping.mapper.label_type_to_id_type[id_type]
                identifier = _identifier

            if id_type == 'mir-pre':

                _identifier = mapping.map_name0(
                    identifier,
                    id_type,
                    'mirbase',
                    ncbi_tax_id = taxon,
                )

                if _identifier and _identifier != identifier:

                    identifier = _identifier
                    id_type = 'mirbase'

        entity_type = entity_type or self._get_entity_type(identifier)

        self.identifier = identifier
        self.id_type = id_type
        self.entity_type = entity_type
        self.taxon = taxon




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


    @classmethod
    def filter_entity_type(cls, entities, entity_type):
        """
        Filters an iterable of entities or identifiers keeping only the ones
        of type(s) in ``entity_type``.
        
        :param iterable entities:
            A list, set, tuple or other iterable yielding entities or
            identifiers.
        :param str,set entity_type:
            One or more entity types e.g. ``{'protein', 'mirna'}``.
        
        :returns:
            Same type of object as ``entities`` if the type of the object is
            list, set or tuple, otherwise a generator.
        """
        
        if not entity_type or not entities:
            
            return entities
        
        entity_type = common.to_set(entity_type)
        obj_type = (
            type(entities)
                if isinstance(entities, common.list_like) else
            lambda x: x
        )
        
        return obj_type(
            e
            for e in entities
            if cls._get_entity_type(e) in entity_type
        )
    
    
    @classmethod
    def only_proteins(cls, entities):
        
        return cls.filter_entity_type(entities, entity_type = 'protein')
    
    
    @classmethod
    def only_complexes(cls, entities):
        
        return cls.filter_entity_type(entities, entity_type = 'complex')
    
    
    @classmethod
    def only_mirnas(cls, entities):
        
        return cls.filter_entity_type(entities, entity_type = 'mirna')


    @classmethod
    def count_entity_type(cls, entities, entity_type):
        """
        Counts elements in an iterable of entities or identifiers of type(s)
        in ``entity_type``.
        
        :param iterable entities:
            A list, set, tuple or other iterable yielding entities or
            identifiers.
        :param str,set entity_type:
            One or more entity types e.g. ``{'protein', 'mirna'}``.
        
        :returns:
            int
        """
        
        entities = (
            entities
                if isinstance(entities, common.list_like) else
            list(entities)
        )
        
        return len(
            cls.filter_entity_type(
                entities,
                entity_type = entity_type,
            )
        )


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


    @classmethod
    def info(cls, identifier):

        if cls._is_protein(identifier):

            return utils_uniprot.info(identifier)


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
