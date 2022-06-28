#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This is a script to exemplarise the usage of an adapter 
(pypath_biocypher_adapter) to create a BioCypher-compatible graph 
database from pypath objects. The mode of interaction with the graph in 
this case is the passing of the locally instantiated Neo4j driver to the 
biocypher module, which resembles just one of multiple possible modes of 
interaction. Merging nodes and edges with an UNWIND APOC call is safe, 
but relatively slow.

Copyright 2021, Heidelberg University Clinic

File author(s): Sebastian Lobentanzer
                ...

Distributed under GPLv3 license, see LICENSE.txt.

In the graph, we have the following components:
- nodes: an id (our primary identifier), a label (:Protein, :Complex 
    ...), a dict of properties
    - NOTE: we also have a property called "label", which is the 
        human-readable id of the node
- edges: source and target id (from primary identifiers), a label
    (:POST_TRANSLATIONAL ...), a dict of properties
- node labels are nouns written in CamelBack, edge labels are verbs 
    written in UPPERCASE, properties are lowercase_with_underscore
    
Todo: 
    * duplicate relationships handling
    * progress bar?
    * batch processing
        - via APOC in parts
        - output csv file for admin import without safety

"""

import pypath.biocypher.adapter as adapter

def write_online():
    # Instantiating the adapter class.
    # We are creating a new database, so we wipe and initialise the 
    # local Neo4j instance. Set `wipe = False` if you want to update an 
    # existing BioCypher graph.
    bcy_adapter = adapter.BiocypherAdapter(wipe=True)

    # Build a pypath network database:
    bcy_adapter.build_python_object()

    # Load the python database object into the connected Neo4j DB:
    bcy_adapter.translate_python_object_to_neo4j()
    # quite slow this one, even only with "network".
    # interactions as nodes is taxing

    # create another adapter without wipe to test meta node 
    # functionality
    # bcy_adapter = adapter.BiocypherAdapter(wipe = False)

def write_offline():
    # Instantiating adapter in offline mode
    bcy_adapter = adapter.BiocypherAdapter(offline=True)

    # Build a pypath network database:
    bcy_adapter.build_python_object()

    # Write CSVs for admin import function
    bcy_adapter.write_to_csv_for_admin_import(db_name="import")
    # TODO translate: unfound entities

def main():

    profile = False
    write_online = False # whether to write to graph via the driver
    write_offline = True # whether to write CSVs for admin import function

    if profile:
        import cProfile, pstats, io
        profile = cProfile.Profile()
        profile.enable()

    if write_online:
        write_online()

    if write_offline:
        write_offline()

    if profile:
        profile.disable()
        s = io.StringIO()
        sortby = pstats.SortKey.CUMULATIVE
        ps = pstats.Stats(profile, stream=s).sort_stats(sortby)
        ps.print_stats()
        # print(s.getvalue())
        filename = "create_network.prof"
        ps.dump_stats(filename)


if __name__ == '__main__':

    main()
