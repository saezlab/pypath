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
BioCypher graph database interface. Loads various pypath objects, such as
databases, resources and records into neo4j, using the BioCypher database
schema.
"""

from __future__ import annotations

import importlib as imp
import os
from typing import Any, Generator, Optional

import biocypher
import neo4j_utils
import yaml

import pypath
import pypath.core.network as network_mod
import pypath.omnipath as op
import pypath.resources.network as netres
import pypath.share.session as _session
from pypath.core import annot
from pypath.core import complex as cmplx

__all__ = ['Adapter']


class Adapter(_session.Logger):
    """
    Channel OmniPath database contents into BioCypher.

    Here we implement a temporary solution for writing OmniPath data into
    BioCypher. It will be replaced by methods integrated into fundamental
    pypath object, especially Resource objects.
    """

    _network_param = {
        'name': 'dummy',
        'module': 'network',
        'args': {
            'resources': (
                netres.pathway['signor'],
                netres.pathway['signalink3'],
                netres.pathway['ca1'],
                netres.transcription['signor'],
                netres.transcription['oreganno'],
                netres.mirna_target['signor'],
                netres.mirna_target['ncrdeath'],
            ),
        },
    }


    def __init__(
            self,
            driver: neo4j_utils.Driver | None = None,
            db_name: str | None = None,
            db_uri: str | None = None,
            db_user: str | None = None,
            db_passwd: str | None = None,
            wipe: bool = False,
            offline: bool = False,
            network: pypath_network.Network | None = None,
            **kwargs
        ):
        """
        Load pypath database obejcts into biocypher (Neo4j).

        Args
            driver:
                A ``neo4j.Driver`` instance, created by, for example,
                ``neo4j.GraphDatabase.driver``.
            db_name:
                Name of the database (Neo4j graph) to use.
            db_uri:
                Protocol, host and port to access the Neo4j server.
            db_user:
                Neo4j user name.
            db_passwd:
                Password of the Neo4j user.
            network:
                A network database object.
            wipe:
                Wipe the database after connection, ensuring the data
                is loaded into an empty database.
            kwargs:
                Passed to ``biocypher.Driver``, and ultimately to
                ``neo4j_utils.Driver``.

        Details:
            The connection can be defined in three ways:
            * Providing a ready ``neo4j.Driver`` instance
            * By URI and authentication data
            * By a YML config file
        """

        _session.Logger.__init__(self, name = 'bcy_adapter')

        self.bcy = biocypher.Driver(
            driver = driver,
            db_name = db_name,
            db_uri = db_uri,
            db_user = db_user,
            db_passwd = db_passwd,
            wipe=wipe,
            offline=offline,
            **kwargs
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
        self.__class__ = new


    def set_network(self, network):

        self._log('Added new network.')
        self.network = network


    def dummy_network(self):
        """
        Obtain a small network just to play with.

        For development, to be removed later.
        """

        op.db.define_dataset(**self._network_param)
        netw = op.db.get_db(*self._network_param['name'])
        self.set_network(netw)


    def translate(
            self,
            network: Optional[pypath_network.Network] = None,
            export_csv: bool = False,
            db_name: Optional[str] = None,
        ) -> Generator[tuple, None, None]:
        """
        Loads a pypath network into the biocypher (Neo4j) backend.

        Args
            network:
                A network database object. If `None`, the value of
                :py:attr:`network` will be used.
            export_csv:
                Write the data into CSV files that can be imported by the
                *neo4j-admin* tool. The default behaviour is to connect to
                the server and insert the data into the database.
            db_name:
                In case of writing CSV files, a CLI call of *neo4j-admin*
                will also be created. The name of the target database is
                part of this call.
        """

        network = network or self.network

        if not network:

            self._log('No network provided.')
            return

        # create id-type tuples for nodes
        # to enable translation between pypath and biocypher notation
        # TODO: other node properties (as dict?)
        def gen_nodes(nodes):

            for n in nodes:

                _id = self._process_id(n.identifier)
                _type = n.entity_type
                props = {'taxon': n.taxon, 'label': n.label}
                yield (_id, _type, props)

        nodes = gen_nodes(network.nodes.values())

        if export_csv:

            self.bcy.write_nodes(nodes, db_name = db_name)

        else:

            self.bcy.add_nodes(nodes = nodes)

        # create id-type tuples for edges
        # to enable translation between pypath and biocypher notation
        # TODO: other edge properties (as dict?)
        def gen_edges(edges):

            for e in edges:

                src = self._process_id(e.id_a)
                tar = self._process_id(e.id_b)
                _type = e.type
                props = {
                    'directed': e.directed,
                    'resources': e.sources,
                    'references': e.references,
                    'effect': e.effect,
                }

                yield (src, tar, _type, props)

        edges = gen_edges(network.generate_df_records(with_references = True))

        if export_csv:

            self.bcy.write_edges(edges, db_name = db_name)
            self.bcy.write_import_call()

        else:

            self.bcy.add_edges(edges = edges)


    def _process_id(self, identifier: str) -> str:
        """
        Replace critical symbols in pypath ids so that neo4j doesn't throw
        a type error.
        """

        return str(identifier).replace('COMPLEX:', 'COMPLEX_')


    def load(self, obj: Any):
        """
        Loads any compatible object into the biocypher (Neo4j) database.

        Args
            obj:
                An object from this module compatible with the current
                adapter. Currently the following database objects are
                supported:
                    * :py:class:`pypath.core.network.Network`
        """

        if hasattr(obj, 'nodes') and hasattr(obj, 'interactions'):

            self.translate_python_object_to_neo4j(network = obj)
