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

import pypath.share.common as common


class Field(object):
    """
    Generic field processor, base class for more specific field processors.
    The default behaviour is to return the value at a specific index from
    each record.
    
    :param int,list idx:
        One or more index or indices.
    :param tuple,list,dict compact:
        Special compact definitions which make it easier to describe field
        definitions in a concise way.
    :param dict mapping:
        Mapping rules for extracting values from the raw field content.
    :param str,NoneType sep:
        Within field separator (to split the field content).
    :param callable method:
        A method for processing - in case the built in ways are not
        satisfying, for greatest flexibility.
    """
    
    def __init__(
            self,
            idx = None,
            compact = None,
            mapping = None,
            sep = None,
            method = None,
            **kwargs
        ):
        
        self.param = locals()
        _ = self.param.pop('self')
        self.param.update(kwargs)
        self.setup()
    
    
    def setup(self):
        
        for k, v in iteritems(param):
            
            setattr(self, k, v)
    
    
    def process(self, record):
        """
        Processes a record and returns the value of the field.
        
        :param list record:
            One record (row) to be processed.
        
        :returns:
            The processed value of the field as a ``FieldContent`` object.
            This object either provides the processed value or iterates
            through values if possibly more than one value available from
            the field.
        """
        
        return FieldContent(record[self.idx])


class FieldContent(object):
    """
    Provides a unified interface for accessing processed field contents
    either as a single value or as an iterable.
    """
    
    
    def __init__(self, content):
        
        self.content = content
    
    
    @property
    def value(self):
        
        return self.content
    
    
    def __iter__(self):
        
        for value in (
            self.content
                if isinstance(self.content, common.list_like) else
            (self.content,)
        ):
            
            yield value
