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
from past.builtins import xrange, range

import imp
import collections

try:
    
    import pybel
    
    if hasattr(pybel, 'ob'):
        
        sys.stdout.write(
            'pypath.bel: You have the `openbabel` module installed '
            'instead of pybel.\n'
            'To be able to use `pybel`, create a virtual env and install '
            'it by `pip install pybel`.\n'
        )
        
        del sys.modules['pybel']
        pybel = None
    
except ModuleNotFoundError:
    
    sys.stdout.write(
        'pypath.bel: module `pybel` not available.\n'
        'You won\'t be able to read or write BEL models.\n'
    )
    
    pybel = None


Relationship = collections.namedtuple(
    'Relationship',
    ('subject', 'predicate', 'object'),
)


class Bel(object):
    """
    Converts pypath objects to BEL format.
    """
    
    def __init__(
            self,
            resource,
        ):
        """
        resource :
            Object to be converted.
            E.g. ``pypath.main.PyPath`` or
            ``pypath.ptm.PtmAggregator`` or
            ``pypath.complex.ComplexAggregator`` or
            ``pypath.network.NetworkResource``.
        """
        
        self.relationships = None
        self.bel = None
        self.resource = resource
    
    
    def reload(self):
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    
    def main(self):
        
        self.resource_to_relationships()
        self.relationships_to_bel()
    
    
    def resource_to_relationships(self):
        
        if hasattr(self.resource, 'graph'):
            # PyPath object
            
            self.resource_to_relationships_graph()
            
        elif hasattr(self.resource, ''):
            # PtmAggregator object
            
            self.resource_to_relationships_enzyme_substrate()
            
        elif hasattr(self.resource, 'complexes'):
            # ComplexAggregator object
            
            self.resource_to_relationships_complex()
            
        elif hasattr(self.resource, 'network'):
            # NetworkResource object
            
            self.resource_to_relationships_network()
    
    
    def resource_to_relationships_graph(self):
        
        pass
    
    
    def resource_to_relationships_enzyme_substrate(self):
        
        pass
    
    
    def resource_to_relationships_complex(self):
        
        raise NotImplementedError
    
    
    def resource_to_relationships_network(self):
        
        raise NotImplementedError
    
    
    def relationships_to_bel(self):
        
        pass
    
    
    def export_relationships(self, fname):
        """
        Exports relationships into 3 columns table.
        
        fname : str
            Filename.
        """
        
        with open(fname, 'w') as fp:
            
            _ = fp.write('Subject\tPredicate\tObject\n')
            
            _ = fp.write(
                '\n'.join(
                    '\t'.join(rel)
                    for rel in self.relationships
                )
            )
    
    
    def export_bel(self, fname):
        """
        Exports the BEL model into file.
        
        fname : str
            Filename.
        """
        
        pass
