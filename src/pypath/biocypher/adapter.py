#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This is a module to exemplarise an adapter class used in creating a 
BioCypher-compatible graph database from pypath objects. It is called from
the "create_and_update_db.py" script in this directory.

Copyright 2021, Heidelberg University Clinic

File author(s): Denés Turei
                Sebastian Lobentanzer
                ...

Distributed under GPLv3 license, see LICENSE.txt.
"""


import os
import yaml
import importlib as imp

import neo4j

import pypath.core.network as pypath_network
import pypath.resources.network as pypath_netres
import pypath.share.session as _session

class ToBioCypher(_session.Logger):

    
    def __init__(
        self,
        db_uri = 'neo4j://localhost:7687',
        db_auth = None,
        config_file = 'db_config.yml',
        network = None,
    ):

        _session.Logger.__init__(self, name = 'pp_tbc')

        self._db_config = {
            'uri': db_uri,
            'auth': db_auth,
        }
        self._config_file = config_file

        self.db_connect()

        if network:

            self.set_network(network)


    def reload(self):
        """
        Reloads the object from the module level.
        """

        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)


    def db_connect(self):

        if not all(self._db_config.values()):

            self.read_config()

        # check for database running?
        self.driver = neo4j.GraphDatabase.driver(**self._db_config)

        self._log('Opened database connection.')


    def read_config(self, section = 'default'):

        if self._config_file and os.path.exists(self._config_file):

            self._log('Reading config from `%s`.' % self._config_file)

            with open(self._config_file, 'r') as fp:

                conf = yaml.safe_load(fp.read())

            self._db_config.update(conf[section])
            self._db_config['auth'] = tuple(self._db_config['auth'])


    def set_network(self, network):

        self._log('Added new network.')
        self.network = network


    def db_close(self):

        self.driver.close()


    def __del__(self):

        self.db_close()


    def load_pypath(self):
        """
        We are loading the pypath network module with the two datasets 
        'pathway' and 'mirna_target' as an example. The dataset is preloaded
        to avoid waiting time when applying multiple database tests.
        """
        n = pypath_network.Network()
        n.load(pypath_netres.pathway)
        n.load(pypath_netres.mirna_target)
        return n
