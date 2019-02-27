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
    
    
    def main(self):
        
        self.resource_to_relationships()
        self.relationships_to_bel()
    
    
    def resource_to_relationships(self):
        
        pass
    
    
    def relationships_to_bel(self):
        
        pass
    
    
    def export(self, fname):
        """
        Exports the BEL model into file.
        
        fname : str
            Filename.
        """
        
        pass
