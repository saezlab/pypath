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

import imp
import collections

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


class Network(session_mod.Logger):
    
    
    def __init__(self, records, hdr = None):
        
        session_mod.Logger.__init__(self, name = 'network')
        
        if isinstance(records, pd.DataFrame):
            
            self.records = records
            
        else:
            
            if not isinstance(records, (list, tuple, pd.np.ndarray)):
                
                records = list(records)
            
            hdr = hdr or records[0]._fields
            
            self.records = pd.DataFrame(records, columns = hdr)
    
    
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
        
        return cls(records = list(pa))
