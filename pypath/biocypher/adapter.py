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
#           Tennur Kılıç
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
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
import pypath.core.network as pypath_network
import pypath.omnipath as op
import pypath.resources.network as pypath_netres
import pypath.share.session as _session
from pypath.core import annot
from pypath.core import complex as cmplx


class BiocypherAdapter(_session.Logger):
    """
    Load pypath database obejcts into biocypher (Neo4j).

    Args:
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


    def __init__(
            self,
            driver: Optional[neo4j_utils.Driver] = None,
            db_name: Optional[str] = None,
            db_uri: Optional[str] = None,
            db_user: Optional[str] = None,
            db_passwd: Optional[str] = None,
            wipe: bool = False,
            offline: bool = False,
            network: Optional[pypath.core.network.Network] = None,
            **kwargs
        ):

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
        setattr(self, '__class__', new)


    def set_network(self, network):

        self._log('Added new network.')
        self.network = network


    # this might be good for development, but should be removed later
    def build_python_object(self):
        """
        Builds a network database with two datasets: 'pathway' and
        'mirna_target'. Intended to be an example. The dataset is preloaded
        to avoid waiting time when applying multiple database tests.
        The resulting network database is stored under the :py:attr:`network`
        attribute.
        """

        ### network
        n = pypath_network.Network()
        # exclude = {'TRIP', 'CellChatDB', 'ncRDeathDB', 'DIP', 'Wojtowicz2020'}
        networks_literal = [
            "pathway", # standard activity flow
            "pathway_noref", # above without literature references
            "interaction", # large undirected databases such as intact
            "ptm", # == "enzyme_substrate" 
            "ptm_noref", # enzyme_substrate without literature references
            "transcription_onebyone", # direct literature curated TF-target interaction
            "dorothea",
            "mirna_target", # direct literature curated miRNA-target interaction
            "tf_mirna", # direct literature curated TF-miRNA interaction
            "lncrna_target", # direct literature curated lncRNA-target interaction
            "ligand_receptor", # all ligand receptor interactions, with and without literature references
            "small_molecule_protein", # all small molecule-protein interactions, with and without literature references 
        ]
        for net in networks_literal:
            n.load(getattr(pypath_netres, net))

            # for each net, build biocypher graph


        ### complexes
        for cls in cmplx.complex_resources:
            # each resource annotates complexes
            complex_db = getattr(complex, cls)()


        ### annotations
        # annot = op.db.get_db("annotations") # get everything at once
        for cls in annot.protein_sources_default:
            # each resource annotates proteins or genes
            annot_db = getattr(annot, cls)()
        
        for cls in annot.complex_sources_default:
            # each resource annotates complexes
            annot_db = getattr(annot, cls)()


        ### cell-cell
        cell_db = op.db.get_db("intercell") # get everything at once
        for role, participants in cell_db.classes.items():
            # import into biocypher for each resource separately
            pass


        ### enzyme-substrate
        enz_db = op.db.get_db("enz_sub") # get everything at once
        for (enz, sub), ptms in enz_db.items():
            # import individual interactions into biocypher 
            pass

        self.set_network(n)


    def translate_python_object_to_neo4j(
            self,
            network: Optional[pypath.core.network.Network] = None,
        ) -> Generator[tuple, None, None]:
        """
        Loads a pypath network into the biocypher (Neo4j) backend.

        Args:
            - network (pypath.core.network.Network): A network database
              object. If `None`, the value of :py:attr:`network` will be
              used.
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
                _props = {"taxon": n.taxon}
                yield (_id, _type, _props)

        id_type_tuples = gen_nodes(network.nodes.values())
        self.bcy.add_nodes(id_type_tuples)

        # create id-type tuples for edges
        # to enable translation between pypath and biocypher notation
        # TODO: other edge properties (as dict?)
        def gen_edges(edges):

            for e in edges:

                _src = self._process_id(e.id_a)
                _tar = self._process_id(e.id_b)
                _type = e.type

                yield (_src, _tar, _type)

        src_tar_type_tuples = gen_edges(network.generate_df_records())
        self.bcy.add_edges(src_tar_type_tuples)


    def write_to_csv_for_admin_import(
            self,
            network: Optional[pypath.core.network.Network] = None,
            db_name: Optional[str] = None,
        ):
        """
        Loads a pypath network into the biocypher (Neo4j) backend using
        the fast Admin Import function, which requires text files that
        need to be properly formatted since it turns off safety measures
        at import.

        Args:
            network:
                A network database object. If `None`, the value of
                :py:attr:`network` will be used.
        """

        network = network or self.network

        if not network:

            self._log("No network provided.")

            return

        # write nodes
        def gen_nodes(nodes):

            for n in nodes:

                _id = self._process_id(n.identifier)
                _type = n.entity_type
                _props = {"taxon": n.taxon, "label": n.label}

                yield (_id, _type, _props)

        id_type_tuples = gen_nodes(network.nodes.values())

        self.bcy.write_nodes(id_type_tuples, db_name=db_name)

        # write edges
        def gen_edges(edges):

            for e in edges:

                _src = self._process_id(e.id_a)
                _tar = self._process_id(e.id_b)
                _type = e.type
                _props = {"effect": e.effect, "directed": e.directed}

                yield (_src, _tar, _type, _props)

        src_tar_type_tuples = gen_edges(network.generate_df_records())

        self.bcy.write_edges(src_tar_type_tuples, db_name=db_name)

        self.bcy.write_import_call()


    def _process_id(self, identifier: str) -> str:
        """
        Replace critical symbols in pypath ids so that neo4j doesn't throw
        a type error.
        """

        return str(identifier).replace('COMPLEX:', 'COMPLEX_')


    def load(self, obj: Any):
        """
        Loads any compatible object into the biocypher (Neo4j) database.

        Args:
            obj:
                An object from this module compatible with the current
                adapter. Currently the following database objects are
                supported:
                    * :py:class:`pypath.core.network.Network`
        """

        if hasattr(obj, 'nodes') and hasattr(obj, 'interactions'):

            self.translate_python_object_to_neo4j(network = obj)
