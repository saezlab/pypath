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

from future.utils import iteritems

import os

try:
    import cPickle as pickle
except:
    import pickle


import pypath.dataio as dataio
import pypath.common as common
import pypath.session_mod as session_mod


class AbstractResource(session_mod.Logger):
    """
    Generic class for downloading, processing and serving
    data from a resource.
    """


    def __init__(
            self,
            name,
            ncbi_tax_id = 9606,
            input_method = None,
            input_args = None,
            dump = None,
            data_attr_name = None,
            **kwargs
        ):
        """
        name : str
            Custom name for the resource.
        input_method : callable
            Method providing the input data.
        """
        
        if not hasattr(self, '_log_name'):
            
            session_mod.Logger.__init__(self, name = 'resource')
        
        self.dump = dump
        self.name = name
        self._data_attr_name = data_attr_name or 'data'
        self._input_method = input_method
        self.input_args = input_args or {}
        self.ncbi_tax_id = ncbi_tax_id


    def load(self):

        self.set_method()
        from_dump = self.from_dump()
        
        if not from_dump:
            
            self.load_data()
            self.process()

        if hasattr(self, 'data'):

            delattr(self, 'data')


    def set_method(self):
        """
        Sets the data input method by looking up in ``dataio`` module if
        necessary.
        """

        if (
            isinstance(self._input_method, common.basestring) and
            hasattr(dataio, self._input_method)
        ):

            self.input_method = getattr(dataio, self._input_method)

        elif callable(self._input_method):

            self.input_method = self._input_method


    def load_data(self):
        """
        Loads the data by calling ``input_method``.
        """
        
        self._log('Loading data from `%s`.' % self.name)
        
        self.set_method()
        
        if hasattr(self, 'input_method'):

            self.data = self.input_method(**self.input_args)


    def process(self):
        """
        Calls the ``_process_method``.
        """
        
        self._log('Processing data from `%s`.' % self.name)
        self._process_method()


    def _process_method(self):

        pass
    
    
    def from_dump(self):
        
        if self.dump is not None:
            
            if (
                isinstance(self.dump, common.basestring) and
                os.path.exists(self.dump)
            ):
                
                with open(self.dump, 'rb') as fp:
                    
                    self._from_dump = pickle.load(fp)
                
            else:
                
                self._from_dump = self.dump
            
            self._from_dump_callback()
            
            return True
        
        return False
    
    
    def _from_dump_callback(self):
        
        if hasattr(self, '_from_dump'):
            
            setattr(self, self._data_attr_name, self._from_dump)
            delattr(self, '_from_dump')
            delattr(self, 'dump')
    
    
    def save_to_pickle(self, pickle_file):
        
        with open(pickle_file, 'wb') as fp:
            
            pickle.dump(
                obj = getattr(self, self._data_attr_name),
                file = fp,
            )


class ResourceAttributes(object):
    
    
    def __init__(
            self,
            name,
            data_type,
            evidence_types = None,
            **kwargs
        ):
        
        self.name = name
        self.data_type = data_type
        self.evidence_types = evidence_types or set()
        
        for attr, value in iteritems(kwargs):
            
            setattr(self, attr, value)
    
    
    def __eq__(self, other):
        
        return (
            self.name == other.name and self.data_type == other.data_type
                if isinstance(other, self.__class__) else
            self.name == other
        )


class NetworkResource(ResourceAttributes):
    
    
    def __init__(
            self,
            name,
            interaction_type = 'PPI',
            data_model = None,
            evidence_types = None,
            **kwargs,
        ):
        
        ResourceAttributes.__init__(
            self,
            name = name,
            data_type = 'network',
            interaction_type = interaction_type,
            evidence_types = evidence_types,
            data_model = data_model,
            **kwargs,
        )
