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
#                  Olga Ivanova
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

import importlib as imp
import collections
import itertools

import pandas as pd

import pypath.session_mod as session_mod
import pypath.interaction as interaction_mod
import pypath.evidence as evidence
import pypath.entity as entity_mod
import pypath.common as common

#Interaction = collections.namedtuple(
    #'Interaction',
    #[
        #'id_a',
        #'id_b',
        #'type_a',
        #'type_b',
        #'directed',
        #'effect',
        #'type',
        #'sources',
        #'references',
    #],
#)
#Interaction.__new__.__defaults__ = (None,) * 7


class Network(session_mod.Logger):
    """
    Represents a molecular interaction network. Provides various methods to
    query the network and its components. Optionally converts the network
    to a ``pandas.DataFrame`` of interactions.
    
    :arg list,dict resources:
        One or more lists or dictionaries containing
        ``pypath.resource.NetworkResource`` objects.
    :arg bool make_df:
        Create a ``pandas.DataFrame`` already when creating the instance.
        If no network data loaded no data frame will be created.
    """
    
    def __init__(
            self,
            resources = None,
            make_df = False,
            df_by_source = False,
            df_with_references = False,
            df_columns = None,
            df_dtype = None,
            **kwargs
        ):
        
        session_mod.Logger.__init__(self, name = 'network')
        
        self._log('Creating network object.')
        
        self.df_by_source = df_by_source
        self.df_with_references = df_with_references
        self.df_columns = df_columns
        self.df_dtype = df_dtype
        
        self.reset()
        self.load(resources = resources, make_df = make_df)
    
    
    def reload(self):
        """
        Reloads the object from the module level.
        """
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    
    def reset(self):
        """
        Removes network data i.e. creates empty interaction and node
        dictionaries.
        """
        
        self.interactions = {}
        self.nodes = {}
        self.nodes_by_label = {}
    
    
    def load(self, resources = None, make_df = False):
        
        pass
    
    
    def make_df(
            self,
            records = None,
            by_source = None,
            with_references = None,
            columns = None,
            dtype = None,
        ):
        """
        Creates a ``pandas.DataFrame`` from the interactions.
        """
        
        self._log('Creating interactions data frame.')
        
        by_source = by_source if by_source is not None else self.df_by_source
        with_references = (
            with_references
                if with_references is not None else
            self.df_with_references
        )
        columns = columns or self.df_columns
        dtype = dtype or self.df_dtype
        
        if not dtype:
            
            dtype = {
                'id_a': 'category',
                'id_b': 'category',
                'type_a': 'category',
                'type_b': 'category',
                'effect': 'int8',
                'type': 'category',
                'sources': 'category' if by_source else 'object',
                'references': 'object' if with_references else 'category',
            },
        
        if not records:
            
            records = self.iter_interactions(
                by_source = by_source,
                with_references = with_references,
            )
        
        if not isinstance(records, (list, tuple, pd.np.ndarray)):
            
            records = list(records)
        
        if not columns and hasattr(records[0], '_fields'):
            
            columns = records[0]._fields
        
        self.df = pd.DataFrame(
            records,
            columns = columns,
            dtype = dtypes,
        )
        
        ### why?
        if dtype:
            
            self.df = self.df.astype(dtype)
        
        self._log(
            'Interaction data frame ready. '
            'Memory usage: %s ' % common.df_memory_usage(self.df)
        )
    
    
    def iter_interactions(self, by_source = False, with_references = False):
        
        pass
    
    
    @classmethod
    def from_igraph(cls, pa, **kwargs):
        """
        Creates an instance from an ``igraph.Graph`` based
        ``pypath.main.PyPath`` object.
        
        :arg pypath.main.PyPath pa:
            A ``pypath.main.PyPath`` object with network data loaded.
        """
        
        obj = cls(**kwargs)
        
        for ia in pa.graph.es['attrs']:
            
            obj.add_interaction(ia)
        
        return obj
    
    
    def add_interaction(self, interaction):
        
        key = (interaction.a, interaction.b)
        
        if key not in interactions:
            
            self.interactions[key] = interaction
            
        else:
            
            self.interactions[key] += interaction
        
        self.add_node(interaction.a)
        self.add_node(interaction.b)
    
    
    def add_node(self, entity):
        
        self.nodes[entity.identifier] = entity
        self.nodes_by_label[entity.label] = entity
    
    
    @property
    def resources(self):
        """
        Returns a set of all resources.
        """
        
        return set.union(*self.df.sources)
    
    
    def entities_by_resource(self):
        """
        Returns a dict of sets with resources as keys and sets of entity IDs
        as values.
        """
        
        return dict(
            (
                resource,
                set(
                    itertools.chain(
                        *self.df[
                            [
                                resource in resources
                                for resources in self.df.sources
                            ]
                        ][['id_a', 'id_b']].values
                    )
                )
            )
            for resource in self.resources
        )
    
    
    def to_igraph(self):
        
        raise NotImplementedError
