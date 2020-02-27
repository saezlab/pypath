#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module.
#  Provides a high level interface for managing builds of the
#  OmniPath databases.
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

import importlib as imp
import os
import itertools
import json

import pypath.share.session as session

_logger = session_mod.Logger(name = 'db_define')
_log = _logger._log
_console = _logger._console


class DatabaseDefinition(object):


    def __init__(self, label, **kwargs):

        self.label = label

        _def = kwargs.pop('def', {}):

        for k, v in itertools.chain(
            iteritems(kwargs),
            iteritems(_def),
        ):

            setattr(self, k, v)


    def __repr__(self):

        return '<%s database: %s>' % (self.class.capitalize(), self.label)


    @classmethod
    def from_json(cls, path, label = None):
        """
        :param str path:
            Path to JSON file with database definition.
        """

        data = cls._parse_json(path = path, label = label)

        return cls(**data)
    
    
    @staticmethod
    def _parse_json(path, label = None):
        
        data = {}
        
        if not os.path.exists(path):
            
            _console('No such file: `%s`.' % path)
            
        else:
            
            data = json.load(path)

            if label:

                if label in data:

                    data = data[label]
                    data['label'] = label

                else:

                    _console(
                        'Entry `%s` not available in file `%s`.' % (
                            label,
                            path,
                        )
                    )
        
        return data


    @classmethod
    def from_dict(cls, dct, label = None):
        """
        :param dict dct:
            Dictionary containing the parameters for the database definition.
        """

        if label:

            dct['label'] = label

        return cls(**dct)


class DatabaseClass(object):
    """
    Describes a class of databases which can be filled with different data
    but here the module and the class implementing the database are defined.
    """

    def __init__(self, label, module, dbclass):

        self.label = label
        self.module = module
        self.dbclass = dbclass
    
    
    def __repr__(self):
        
        return (
            '<Database class `%s`, module: `%s`, class or method: `%s`>' % (
                self.label,
                self.module,
                self.dbclass.__name__
                    if hasattr(self.dbclass, '__name__') else
                self.dbclass,
            )
        )


    def get_class(self):

        if callable(self.dbclass):

            return self.dbclass

        else:

            try:

                mod = importlib.import_module(self.module)
                
                if hasattr(mod, self.dbclass):
                    
                    return getattr(mod, self.dbclass)
                    
                else:
                    
                    _console(
                        'Module `%s` has no class or method `%s`.' % (
                            self.module,
                            self.dbclass,
                        )
                    )
                
            except ImportError:
                
                _console('Failed to import `%s`.' % self.module)
    
    
    @classmethod
    def from_json(self, path, label = None):
        
        data = DatabaseDefinition._parse_json(path = path, label = label)
        
        return cls(**data)
    
    
    @classmethod
    def from_dict(cls, dct):
        
        return cls(**dct)
