#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This is a module to exemplarise an adapter class used in creating a 
BioCypher-compatible graph database from pypath objects. It is called from
the "create_and_update_db.py" script in this directory.

Copyright 2021, Heidelberg University Clinic

File author(s): Dénes Türei
                Sebastian Lobentanzer
                ...

Distributed under GPLv3 license, see LICENSE.txt.
"""


import os
import yaml
import importlib as imp

import biocypher

import pypath.core.network as pypath_network
import pypath.resources.network as pypath_netres
import pypath.share.session as _session


class BiocypherAdapter(_session.Logger):


    def __init__(
        self,
        driver = None,
        db_uri = 'neo4j://localhost:7687',
        db_auth = None,
        config_file = 'db_config.yml',
        network = None,
    ):

        _session.Logger.__init__(self, name = 'bcy_adapter')

        self.bcy = biocypher.Driver(
            driver = driver,
            db_uri = db_uri,
            db_auth = db_auth,
            config_file = config_file,
        )

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


    def set_network(self, network):

        self._log('Added new network.')
        self.network = network


    def build_network(self):
        """
        Builds a network database with two datasets: 'pathway' and
        'mirna_target'. Intended to be an example. The dataset is preloaded
        to avoid waiting time when applying multiple database tests.
        The resulted network database is stored under the :py:attr:`network`
        attribute.
        """

        n = pypath_network.Network()
        n.load(pypath_netres.pathway)
        n.load(pypath_netres.mirna_target)

        self.set_network(n)


    def load_network(self, network = None):
        """
        Loads a network into the biocypher (Neo4j) backend.

        Args:
            network (pypath.core.network.Network): A network database object.
                If `None`, the value of :py:attr:`network` will be used.
        """

        network = network or self.network

        if not network:

            self._log('No network provided.')
            return

        self.bcy.add_nodes(network.nodes.values())
        self.bcy.add_edges(network.generate_df_records())
