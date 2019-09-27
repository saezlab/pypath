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

import imp
import collections
import itertools

import pandas as pd

import pypath.session_mod as session_mod


Interaction = collections.namedtuple(
    'Interaction',
    [
        'id_a',
        'id_b',
        'type_a',
        'type_b',
        'directed',
        'effect',
        'type',
        'sources',
        'references',
    ],
)
Interaction.__new__.__defaults__ = (None,) * 7


class Network(session_mod.Logger):
    
    
    def __init__(self, records, dtypes = None, **kwargs):
        
        session_mod.Logger.__init__(self, name = 'network')
        
        if isinstance(records, pd.DataFrame):
            
            self.records = records
            
        else:
            
            if not isinstance(records, (list, tuple, pd.np.ndarray)):
                
                records = list(records)
            
            if 'columns' not in kwargs and hasattr(records[0], '_fields'):
                
                kwargs['columns'] = records[0]._fields
            
            self.records = pd.DataFrame(records, **kwargs)
        
        if dtypes:
            
            self.records = self.records.astype(dtypes)
    
    
    def reload(self):
        """
        Reloads the object from the module level.
        """
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    
    @classmethod
    def from_igraph(cls, pa):
        """
        Creates an instance using a ``pypath.main.PyPath`` object.
        """
        
        return cls(
            records = pa.__iter__(),
            dtypes = {
                'id_a': 'category',
                'id_b': 'category',
                'type_a': 'category',
                'type_b': 'category',
                'effect': 'int8',
            },
        )
    
    
    @property
    def resources(self):
        """
        Returns a set of all resources.
        """
        
        return set.union(*self.records.sources)
    
    
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
                        *self.records[
                            [
                                resource in resources
                                for resources in self.records.sources
                            ]
                        ][['id_a', 'id_b']].values
                    )
                )
            )
            for resource in self.resources
        )
