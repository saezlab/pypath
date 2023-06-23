#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright 2014-2023
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: see the file `README.rst`
#  Contact: Dénes Türei (turei.denes@gmail.com)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      https://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: https://pypath.omnipathdb.org/
#

"""
Classes for extracting values from network specific fields in a record (row)
of raw or preprocessed data.
"""

import pypath.reader.field as field


class DirectionBase(field.Field):
    
    
    def __init__(self, *args, **kwargs):
        
        field.Field.__init__(self, *args, **kwargs)



class Direction(DirectionBase):
    """
    Processes the fields describing the direction of an interaction in a
    record describing a network interaction.
    """
    
    
    def __init__(self, *args, **kwargs):
        
        DirectionBase.__init__(self, *args, **kwargs)


class Sign(DirectionBase):
    """
    Processes the fields describing the effect sign of an interaction in a
    record describing a network interaction.
    """
    
    
    def __init__(self, *args, **kwargs):
        
        DirectionBase.__init__(self, *args, **kwargs)


class Filter(field.Field):
    """
    Processes a field into a positive or negative filter.
    """
    
    def process(self, record):
        
        field.Field.process(self, record)
