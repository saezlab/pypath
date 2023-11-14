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

from future.utils import iteritems

import importlib
import os
import itertools
import json

import pypath.share.session as session
import pypath.omnipath.databases.build as build

_logger = session.Logger(name = 'db_define')
_log = _logger._log
_console = _logger._console


class DatabaseDefinition(object):


    def __init__(self, label, **kwargs):

        self.label = label

        _def = kwargs.pop('def', {})

        for k, v in itertools.chain(
            iteritems(kwargs),
            iteritems(_def),
        ):

            setattr(self, k, v)


    def __repr__(self):

        return '<%s database: %s>' % (self.dbclass.capitalize(), self.label)


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

        data = DatabaseDefinition._read_json(path) or {}

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


    @staticmethod
    def _read_json(path):

        if not os.path.exists(path):

            _console('No such file: `%s`.' % path)

        else:

            with open(path) as json_file:

                return json.load(json_file)


    @classmethod
    def from_dict(cls, dct, label = None):
        """
        :param dict dct:
            Dictionary containing the parameters for the database definition.
        """

        if label:

            dct['label'] = label

        return cls(**dct)


    def get(self, attr):

        if hasattr(self, attr):

            return getattr(self, attr)


class DatabaseClass(object):
    """
    Describes a class of databases which can be filled with different data
    but here the module and the class implementing the database are defined.
    """

    def __init__(self, module, method, label = None):

        self.label = label
        self.module = module
        self.method = method


    def __repr__(self):

        return (
            '<Database class `%s`, module: `%s`, class or method: `%s`>' % (
                self.label,
                self.module,
                self.method.__name__
                    if hasattr(self.method, '__name__') else
                self.method,
            )
        )


    def get_class(self):

        if callable(self.method):

            return self.method

        else:

            try:

                mod = importlib.import_module(self.module)

                if hasattr(mod, self.method):

                    return getattr(mod, self.method)

                else:

                    _console(
                        'Module `%s` has no class or method `%s`.' % (
                            self.module,
                            self.method,
                        )
                    )

            except ImportError:

                _console('Failed to import `%s`.' % self.module)


    @classmethod
    def from_json(self, path, label = None):

        data = DatabaseDefinition._parse_json(path = path, label = label)

        return cls(**data)


    @classmethod
    def from_dict(cls, dct, label = None):

        if label:

            dct['label'] = label

        return cls(**dct)


class DatabaseDefinitionManager(session.Logger):


    def __init__(self, classes = None, databases = None):

        session.Logger.__init__(self, name = 'db_define')

        self._classes = classes or self._default_json('classes')
        self._databases = databases or self._default_json('builtins')

        self.load()


    def __repr__(self):

        return '<Database definitions: %u classes and %u definitions>' % (
            len(self.classes),
            len(self.databases),
        )


    def load(self):

        if isinstance(self._classes, str):

            self._log('Reading database classes from `%s`' % self._classes)
            self._classes = DatabaseDefinition._read_json(self._classes)

        if isinstance(self._databases, str):

            self._log(
                'Reading database definitions from `%s`' % self._databases
            )
            self._databases = DatabaseDefinition._read_json(self._databases)

        self.classes = dict(
            (
                label,
                DatabaseClass.from_dict(label = label, dct = param)
            )
            for label, param in iteritems(self._classes)
        )

        self.databases = dict(
            (
                label,
                DatabaseDefinition.from_dict(label = label, dct = param)
            )
            for label, param in iteritems(self._databases)
        )


    def get_db_class(self, label):

        if label in self.classes:

            return self.classes[label]

        else:

            self._log('No such database class: `%s`.' % label)


    def get_db_definition(self, label):

        if label not in self.databases:

            self._log(
                'Warning: no parameters for label `%s`, '
                'returning empty dict.' % label
            )

        return self.databases[label] if label in self.databases else {}


    def get_class(self, label):

        dbclass = self.get_db_class(label)

        if dbclass:

            return dbclass.get_class()


    def class_and_param(self, label):
        """
        For a database definition label returns the class or method and its
        arguments which are necessary to build the database according to
        the definition.
        """

        db_def = self.get_db_definition(label)

        if db_def:

            db_class = db_def.dbclass

            if not callable(db_class):

                if isinstance(db_class, dict):

                    db_class = DatabaseClass(**db_class)

                elif isinstance(db_class, str):

                    db_class = self.get_db_class(db_class)

        return db_class, db_def


    def build(self, label):
        """
        For a database definition label returns an instance of the database:
        creates an instance of the class or calls the method with the
        arguments in the database definition. Returns the database instance.
        """

        db_class, db_def = self.class_and_param(label)

        if db_class:

            return build.build(db_class, db_def)


    @staticmethod
    def _default_json(name):

        return os.path.join(
            session.session().module_root,
            'omnipath',
            'databases',
            '%s.json' % name,
        )
