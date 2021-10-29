#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This is a script to exemplarise the usage of an adapter 
(pypath_biocypher_adapter) to create a BioCypher-compatible graph database 
from pypath objects. The mode of interaction with the graph in this case is
the passing of the locally instantiated Neo4j driver to the biocypher module, 
which resembles just one of multiple possible modes of interaction. Merging
nodes and edges with an UNWIND APOC call is safe, but relatively slow.

Copyright 2021, Heidelberg University Clinic

File author(s): Sebastian Lobentanzer
                ...

Distributed under GPLv3 license, see LICENSE.txt.

In the graph, we have the following components:
- nodes: an id (our primary identifier), a label (:Protein, :Complex ...), a
    dict of properties
    - NOTE: we also have a property called "label", which is the human-readable
        id of the node
- edges: source and target id (from primary identifiers), a label
    (:POST_TRANSLATIONAL ...), a dict of properties
- node labels are nouns written in CamelBack, edge labels are verbs written in
    UPPERCASE, properties are lowercase_with_underscore
    
Todo: 
    * duplicate relationships handling
    * progress bar?
    * batch processing
        - via APOC in parts
        - output csv file for admin import without safety

"""

# Imports
# =======
import pypath.biocypher.adapter as adapter

# BioCypher import
# The BioCypher prototype needs to be installed locally, please follow the 
# instructions at https://github.com/saezlab/BioCypher
import biocypher.driver as bcy_driver


def main():
    # Setup
    # =====

    # instantiate adapter class
    # loads pypath;
    py_tb = adapter.ToBioCypher()

    # hand over driver to biocypher
    db = bcy_driver.DatabaseToNeo4j(py_tb.driver)

    # We are creating a new database, so we wipe and initialise the local Neo4j
    # instance. Skip this if you want to update an existing BioCypher graph.
    db.init_db()


    # Data input
    # ==========

    # loading pypath network object, to avoid waiting for multiple tests (?)
    # may take a while if it is the first time due to downloads
    net = py_tb.load_pypath()

    # create nodes
    db.add_nodes_to_graph(net.nodes.values())

    # create edges
    db.add_edges_to_graph(net.generate_df_records())


if __name__ == '__main__':

    main()
