#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2020
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

"""
Classes for extracting values from fields in a record (row) of raw or
preprocessed data.
"""

from future.utils import iteritems


class Field(object):
    """
    Generic field processor, base class for more specific field processors.
    The default behaviour is to return the value at a specific index from
    each record.
    """
    
    def __init__(
            self,
            idx,
            mapping = None,
        ):
        
        print(locals())
        self.setup(locals())
    
    
    def setup(self, param):
        
        for k, v in iteritems(param):
            
            setattr(self, k, v)
    
    
    def process(self, record):
        """
        Processes a record and returns the value of the field.
        
        :param list record:
            One record (row) to be processed.
        
        :returns:
            The processed value of the field.
        """
        
        return record[self.idx]
