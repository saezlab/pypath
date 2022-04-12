#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2022
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: Dénes Türei (turei.denes@gmail.com)
#           Nicolàs Palacio
#           Sebastian Lobentanzer
#           Erva Ulusoy
#           Olga Ivanova
#           Ahmet Rifaioglu
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

"""
This file is part of the `pypath` python module. It provides user
settings to the remainder of PyPath modules. Settings are gathered from
`pypath/data/settings.yaml`.
"""

from future.utils import iteritems

import os, yaml
import collections
import contextlib

__all__ = ['Settings', 'settings', 'get', 'setup', 'context']

ROOT = os.path.join(
    *os.path.split(
        os.path.abspath(
            os.path.dirname(__file__)
        )
    )[:-1]
)

# import settings from yaml
# we are importing from module data
_settings_yaml = os.path.join(ROOT, 'data', 'settings.yaml')


with open(_settings_yaml, 'r') as f:
    # log? ## logger is not available here, as logging parameters
    # are defined in the settings
    # TODO: would be better to use an immutable dict-like type?
    _from_config_file = yaml.load(f, Loader = yaml.FullLoader)


class Settings(object):
    """
    Class to provide settings to other modules.

    Args:
        **kwargs: key-value pairs to be included in the settings dict
    """

    # TODO put these custom collections into YAML as well?
    # yes, and most of these are deprecated, these were module
    # data files that we got rid of long time ago
    in_datadir = {
        'acsn_names',
        'alzpw_ppi',
        'goose_annot_sql',
        'webpage_main',
        'nrf2ome',
        'ppoint',
        'slk3_nodes',
        'acsn',
        'arn',
        'goose_ancest_sql',
        'goose_terms_sql',
        'lmpid',
        'nci_pid',
        'old_dbptm',
        'slk3_edges',
        'slk01human',
        'deathdomain',
        'license_dir',
    }

    in_cachedir = {
        'pubmed_cache',
        'trip_preprocessed',
        'hpmr_preprocessed',
    }

    in_secrets_dir = {
        'license_secret',
    }

    # special directories with built in default at user level
    pypath_dirs = (
        ('cachedir', 'cache'),
        ('pickle_dir', 'pickles'),
        ('secrets_dir', 'secrets'),
    )

    def __init__(self, _dict = None, **kwargs):

        self.reset_all()
        self.setup(_dict, **kwargs)


    def reset_all(self):
        """
        Main method of updating the settings object from data directory
        structure and the YAML file contents.
        """

        self._settings = {}
        self._context_settings = []


        self.setup(
            dict(
                (
                    k,
                    os.path.join(ROOT, 'data', val)
                        if k in self.in_datadir else
                    val
                )
                for k, val in _from_config_file.items()
            )
        )

        # runtime attributes
        # base directory
        self.setup(basedir = ROOT)

        self.setup(
            dict(
                (
                    _key,
                    os.path.join(
                        os.path.expanduser('~'),
                        '.pypath',
                        _dir,
                    )
                )
                for _key, _dir in self.pypath_dirs
                if self.get(_key) is None
            )
        )

        self.setup(
            dict(
                (
                    k,
                    os.path.join(self.get('cachedir'), _from_config_file[k])
                )
                for k in self.in_cachedir
            )
        )

        self.setup(
            dict(
                (
                    k,
                    os.path.join(self.get('secrets_dir'), _from_config_file[k])
                )
                for k in self.in_secrets_dir
            )
        )


    def setup(self, _dict = None, **kwargs):
        """
        Set the values of various parameters in the settings.

        Args:
            _dict: A `dict` of parameters, keys are the option names, values
                are the values to be set.
            kwargs: Alternative way to provide parameters, argument names
                are the option names, values are the corresponding values.

        Returns:
            None
        """

        _dict = self._dict_and_kwargs(_dict, kwargs)
        self._settings.update(_dict)


    @staticmethod
    def _dict_and_kwargs(_dict, kwargs):

        _dict = _dict or {}
        _dict.update(kwargs)

        return _dict


    def get(self, param, override = None, default = None):
        """
        Retrieves the current value of a parameter.

        :param str param:
            The key for the parameter.
        :param object,NoneType override:
            If this value is not None it will be returned instead of the
            settings value. It is useful if the parameter provided at the
            class or method level should override the one in settings.
        :param object,NoneType default:
            If no value is available for the parameter in the current
            settings, this default value will be returned instead.
        """

        return (
            (self[param] if param in self else default)
                if override is None else
            override
        )


    def get_default(self, param):
        """
        Returns the default value of a parameter. If no default value is
        defined for the parameter, None will be returned.

        Args:
            param (str): The name of the parameter.

        Returns:
            The default value of the parameter or None.
        """

        return _defaults.get(param)


    @property
    def contexts(self):

        return reversed(self._context_settings)


    def _from_context(self, param):

        for ctx in self.contexts:

            if param in ctx:

                return ctx[param]


    def _in_context(self, param):

        return any(
            param in ctx
            for ctx in self.contexts
        )


    @contextlib.contextmanager
    def context(self, _dict = None, **kwargs):
        """
        Temporarily alter the values of certain parameters. At exiting the
        context, the original values will be restored. Multiple contexts
        can be nested within each other.

        Args:
            _dict: A `dict` of parameters, keys are the option names, values
                are their values.
            kwargs: Alternative way to provide parameters, argument names
                are the option names, values are the corresponding values.
        """

        try:

            ctx = self._dict_and_kwargs(_dict, kwargs)
            self._context_settings.append(ctx)
            yield

        finally:

            self._context_settings = self._context_settings[:-1]


    @property
    def _numof_contexts(self):

        return len(self._context_settings)


    @property
    def _innermost_context(self):

        return self._context_settings[-1] if self._context_settings else None


    def reset(self, param):
        """
        Reset the parameters to their default values.

        Args:
            param: The name of the parameter to be reset.

        Returns:
            None
        """

        self.setup({param: self.get_default(param)})


    def __getattr__(self, attr):

        if attr in self:

            return self[attr]

        elif attr in self.__dict__:

            return self.__dict__[attr]

        else:

            raise AttributeError(
                '\'%s\' object has no attribute \'%s\'' % (
                    self.__class__.__name__,
                    str(attr)
                )
            )


    def __dir__(self):

        keys = object.__dir__(self)
        [keys.extend(ctx.keys()) for ctx in self.contexts]
        keys.extend(self._settings.keys())
        keys = sorted(set(keys))

        return keys


    def __contains__(self, param):

        return (
            self._in_context(param) or
            param in self._settings
        )


    def __getitem__(self, key):

        if self._in_context(key):

            return self._from_context(key)

        elif key in self._settings:

            return self._settings[key]

        else:

            return None


    def __setitem__(self, key, value):

        self._settings[key] = value


_defaults = Settings()
settings = Settings()


def get(param, override = None):
    """
    The current value of a settings parameter.

    Args:
        param (str): Name of a parameter.
        override: Override the currently valid settings,
            return this value instead.
        default: If no value is set up for the key requested,
            use this default value instead.

    Wrapper of `Settings.get()`.
    """
    return settings.get(param, override = override)


def setup(_dict = None, **kwargs):
    """
    Set the values of various parameters in the settings.

    Args:
        _dict: A `dict` of parameters, keys are the option names, values
            are the values to be set.
        kwargs: Alternative way to provide parameters, argument names
            are the option names, values are the corresponding values.

    Returns:
        None
    """

    return settings.setup(_dict, **kwargs)


def context(_dict = None, **kwargs):

    return settings.context(_dict, **kwargs)
