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
#           Olga Ivanova
#           Sebastian Lobentanzer
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

ROOT = os.path.join(
    *os.path.split(
        os.path.abspath(
            os.path.dirname(__file__)
        )
    )[:-1]
)

# import settings from yaml
# we are importing from module data
settings_yaml = os.path.join(ROOT, 'data', 'settings.yaml')


with open(settings_yaml, 'r') as f:
    # log? ## logger is not available here, as logging parameters
    # are defined in the settings
    # TODO: would be better to use an immutable dict-like type?
    _defaults = yaml.load(f, Loader = yaml.FullLoader)


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

    def __init__(self, **kwargs):

        self.__dict__.update(kwargs)
        self.reset_all()


    def reset_all(self):
        """
        Main method of updating the settings object from data directory
        structure and the YAML file contents.
        """

        for k, val in _defaults.items():

            if k in self.in_datadir:

                val = os.path.join(ROOT, 'data', val)

            setattr(self, k, val)

        # runtime attributes
        # base directory
        setattr(self, 'basedir', ROOT)

        for _key, _dir in self.pypath_dirs:

            if getattr(self, _key) is None:

                setattr(
                    self,
                    _key,
                    os.path.join(
                        os.path.expanduser('~'),
                        '.pypath',
                        _dir,
                    )
                )

        for k in self.in_cachedir:

            setattr(self, k, os.path.join(self.cachedir, _defaults[k]))

        for k in self.in_secrets_dir:

            setattr(self, k, os.path.join(self.secrets_dir, _defaults[k]))


    def setup(self, **kwargs):
        """
        This function takes a dictionary of parameters and values and sets them
        as attributes of the settings object.

        Args:
        **kwargs: key-value pairs to set in the `settings` object

        Returns:
        None
        """

        for param, value in iteritems(kwargs):

            setattr(self, param, value)


    def get(self, param, value = None):
        """
        Retrieves the current value of a parameter.

        :param str param:
            The key for the parameter.
        :param object,NoneType value:
            If this value is not None it will be returned instead of the settings
            value. It is useful if the parameter provided at the class or method
            level should override the one in settings.
        """

        if value is not None:

            return value

        if hasattr(self, param):

            return getattr(self, param)



    def get_default(self, param):
        """
        Returns the value of the parameter in the defaults object if it
        exists, otherwise returns None.

        Args:
        param: keyword to look for in `defaults`

        Returns:
        The value of the parameter or None.
        """

        if hasattr(defaults, param):

            return getattr(defaults, param)


    def reset(self, param):
        """
        Reset the parameters to their default values.

        Args:
        param: the name of the parameter to be set

        Returns:
        None
        """

        self.setup(**{param: self.get_default(param)})


settings = Settings()

def get(param, value = None):
    """
    Wrapper of Settings.get().
    """
    return settings.get(param, value)


def setup(**kwargs):

    return settings.setup(**kwargs)
