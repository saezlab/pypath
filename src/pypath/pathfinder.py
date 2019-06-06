#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2019
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#                  Nicolàs Palacio
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

"""
PyPath finder

Tool based on OmniPath/PyPath in order to build prior knowledge networks (PKN)
from a list of nodes set by the user.

Copyright (C) 2017 Nicolàs Palacio

Contact: palacio@combine.rwth-aachen.de

Licensed under GNU-GLPv3:
This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

A full copy of the GNU General Public License can be found on:
http://www.gnu.org/licenses/.
"""

# SOME IDEAS:
# Make coffee
# Plots, graphs and stuff.
# Create a results folder to put all the stuff?
# More prints? (that nobody will read because verbose is for debuggers)
# And/or multi-level verbose mode?

#---------------------------------- MODULES ---------------------------------#
from __future__ import print_function
from future.utils import iteritems
from past.builtins import xrange, range, reduce

import sys
import os
import time
import igraph
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

class silent(object):
   
    def __init__(self):
       
        pass
   
    def __enter__(self):
        self.aux = sys.stdout
        sys.stdout = open(os.devnull, 'w')
   
    def __exit__(self, exception_type, exception_value, traceback):
       
        sys.stdout.close()
        sys.stdout = self.aux

# Silent import of PyPath
try:
    with silent():
       
        import pypath
   
except ImportError as err:
    raise err

#------------------------------- INITIALIZATION -------------------------------#
__all__ = ['Main', 'Main2', 'convert_names', 'read_sif', 'sys', 'os', 'time',
           'pypath', 'igraph', 'np', 'pd', 'plt']

#---------------------------------- CLASSES -----------------------------------#
class Main(object): # XXX: DEPRECATE?
    '''
    Creates/loads a PyPath instance (see the package reference for further
    details). By default edges with only high-throughput references are
    removed. Also loads the gene symbols mapping to UniProt IDs, pathway
    attributes and edge directions. If the instance is not loaded from a
    file (filename not provided), user is prompted to input a new file
    name/path to save the pickle (if none is provided, it does nothing).
    Another option is to use an already existing PyPath instance by passing it
    with the keyword argument pa.

    * Attributes:
        - pa [pypath.main.PyPath]: An initialized PyPath instance (see the
          package reference for further details).
        - mappr [pypath.mapping.Mapper]: PyPath mapping object. Helps
          converting the IDs between gene names and UniProt IDs and vice versa.
        - inputs [list]: Contains the names [str] of the nodes that are
          considered as inputs. This means that nodes contained in this list
          will be considered network inputs and will not have incoming edges in
          the final PKN.
        - query [dict]: The node query containing all the nodes that will form
          the network. Key/value pairs consist of protein names [str], and
          their corresponding UniProt ID [str].
        - nodes [pandas.DataFrame]: Contains some information about the nodes
          in the network. More specifically: Name, UniProtID, Indegree,
          Outdegree and Input (True if it is, False otherwise).
        - edges [pandas.DataFrame]: Contains some information about the edges
          in the network, like: Pair (source-target nodes as indices), Name
          (like previous but as node names), Refs and Type (of relation, 1 for
          stimulation and -1 for inhibition). Please note that the last
          information should be confirmed in the literature manually as it can
          be ambiguous, incomplete or incorrect (See Interaction (+, -) column).
        - network [igraph.Graph]: Contains the PKN as an igraph object. See the
          corresponding documentation for more information about object methods
          and attributes.
        - mode [str]: Specifies the type of identifiers used on the input.
          Options are 'genesymbol' for the gene names and 'uniprot' for UniProt
          IDs.
    '''

    def __init__(self, filename=None, pa=None):
        '''
        This function creates/modifies the following attributes:
            - pa, mappr, query.

        * Arguments:
            - filename [str]: Optional, None by default. File name or path to
              the file where the PyPath pickle should be loaded. If none is
              provided, it loads all the data from scratch.
            - pa [pypath.main.PyPath]: Optional, None by default. An
              initialized PyPath instance (see the package reference for
              further details).
        '''

        self.query = {}
        self.mappr = pypath.mapping.Mapper(9606)

        if pa:
            self.pa = pa

        else:
            self.pa = pypath.PyPath()

            if filename:
                self.pa.init_network(pfile=filename)

            else:
                self.pa.init_network()

            self.pa.remove_htp()
            self.pa.genesymbol_labels()
            self.pa.pathway_attributes()
            self.pa.load_all_pathways()
            self.pa.get_directed()

            if not filename:
                filename = raw_input('Enter name/path to save the pickle: ')
                if filename is not '':
                    self.pa.save_network(filename)
                    print('PyPath instance saved at %s' %filename)

    def build_network(self, m_nodes, i_nodes, max_in=5, max_out=5, breakif=5,
                      verbose=False, mode='genesymbol'):
        '''
        Builds up the PKN from a set of nodes. We discriminate between the
        input and measured nodes. The former ones do not need to have incoming
        edges, the latter contains the rest of nodes that need to be in the PKN.

        This function creates/modifies the following attributes:
            - query, nodes, edges, network, inputs, mode.

        * Arguments:
            - m_nodes [list]: List of measured/other nodes that need to be in
              the PKN and are not inputs. The names [str] in the list should
              match UniProt names for everything to run automatycally.
              Otherwise, the program will ask for manual input of the
              corresponding UniProt ID.
            - i_nodes [list]: Similar to m_nodes, but the nodes included in
              this list will be considered inputs and will not have incoming
              edges in the final PKN.
            - max_in [int]: Optional, 5 by default. Determines the maximum
              number of incoming edges for a node in order to reduce the
              network complexity.
            - max_out [int]: Optional, 5 by default. Determines the maximum
              number of outgoing edges for a node in order to reduce network
              complexity.
            - breakif [int]: Optional, 5 by default. Sets a threshold for which
              edges whose number of references are equal or above, they cannot
              be removed.
            - verbose [bool]: Optional, False by default. Sets on/off the
              printing of messages.
            - mode [str]: Optional, 'genesymbol' by default. Specifies the type
              of identifiers used on the input. Options are 'genesymbol' for
              the gene names and 'uniprot' for UniProt IDs.
        '''

        nodes = set(list(m_nodes) + list(i_nodes))
        self.inputs = i_nodes
        self.mode = mode

        if verbose:
            print('Initial requested query has %d elements.' %len(nodes))
            print('Expanding query...')

        self.build_query(nodes, verbose=verbose)

        if verbose:
            print('Expanding query...[DONE]\n')

        self.network = self.pa.dgraph.induced_subgraph(self.query.values(),
                                           implementation='create_from_scratch')

# FIXME: Some nodes still get disconnected from the network despite having
#        direct neighbors in the query

        self.update_tables()

        if verbose:
            print('Former network built:')
            print('{} nodes and {} edges.'.format(self.network.vcount(),
                                                  self.network.ecount()))
            print('Reducing network...')

# ----->
# XXX: remove after debug
        self.pre_reduction = self.network.copy()
        #igraph.plot(self.pre_reduction, 'pre_reduction.pdf')
# <-----

        self.reduce_inout(max_in=max_in, max_out=max_out, breakif= breakif,
                          verbose=verbose)

        if verbose:
            print('Reducing network...[DONE]\n')
            print('New network:',)
            print('{} nodes and {} edges.'.format(self.network.vcount(),
                                                  self.network.ecount()))

    def build_query(self, nodes, verbose=False):
        '''
        For each node, checks if there exists a direct connection to any other
        node in the set. If no direct connection is found, looks for them in
        the nearby neighbourhood and adds all the intermediate nodes that make
        these connections. If this is not the case, further neighbours are
        checked until at least one path is found.

        This function modifies the following attributes:
            - query.

        * Arguments:
            - nodes [list]: List of nodes that need to be in the PKN. The names
              [str] in the list should match UniProt names for everything to
              run automatycally. Otherwise, the program will ask for manual
              input of the corresponding UniProt ID.
            - verbose [bool]: Optional, False by default. Sets on/off the
              printing of messages.
        '''

        self.update_query(nodes)

        to_add = set()

        # Add all intermediate nodes so that all nodes in query are
        # connected

        for name, n in iteritems(self.query):
            if name in self.inputs or n in self.inputs:
                neighbors = self.pa.dgraph.vs(self.pa.dgraph.neighbors(n,
                                              mode='OUT'))['name']

            else:
                neighbors = self.pa.dgraph.vs(self.pa.dgraph.neighbors(n,
                                              mode='ALL'))['name']

            direct = [i for i in neighbors if i in self.query.itervalues()]

            if direct:
                if verbose:
                    print('Node {} has {} connections in query'.format(name,
                                                                    len(direct)))

            else:
                if verbose:
                    print('No direct connections found for', name)

                paths = []

                for mname, m in iteritems(self.query):
                    if m == n:
                        pass

                    else:
                        if mname not in self.inputs and m not in self.inputs:
                            paths.append(self.pathfinder(n, m))

                        if name not in self.inputs and n not in self.inputs:
                            paths.append(self.pathfinder(m, n))
                       
                lens = np.array([len(i) for i in paths])

                aux = []
                top = 1
                pending = True
# ----->
# XXX: remove after debug
                #print(top)
                #print(lens)
# <-----

                if sum(lens) == 0:
                    pending = False
                    print('\t***WARNING: No connections found for %s***' % name)

                while pending:
                    if top in lens:
                        if verbose:
                            print('\tFound %d shortest paths to query' % np.sum(
                                                                   lens == top))
                            print('with %d intermediates' %top)

                        aux = np.unique([paths[i] for i in
                                         np.where(lens == top)[0]])
                        pending = False
                       
                    else:
                        top += 1
                       
                to_add.update(aux)

        # Check that all nodes that are not an input have at least one incoming
        # edge
        if verbose:
            print('\nChecking that all non-input nodes have at least one',)
            print('incoming edge:')

        temp_network = self.pa.dgraph.induced_subgraph(self.query.values(),
                                           implementation='create_from_scratch')
        temp_nodes = pd.DataFrame(index=range(temp_network.vcount()),
                                  columns=['Name', 'UniProtID', 'Indegree',
                                           'Input'])
        temp_nodes['Name'] = temp_network.vs['label']
        temp_nodes['UniProtID'] = temp_network.vs['name']
        temp_nodes['Indegree'] = temp_network.indegree()
        if self.mode == 'genesymbol':
            temp_nodes['Input'] = [i in self.inputs for i in temp_nodes['Name']]

        elif self.mode == 'uniprot':
            temp_nodes['Input'] = [i in self.inputs for i in
                                   temp_nodes['UniProtID']]

        else:
            raise Exception('ERROR: Invalid mode selection. ' +
                            'Please use "genesymbol" or "uniprot"')

        for row in temp_nodes[temp_nodes.loc[:, 'Indegree'] == 0].values:
            name, n = row[:-2]
            if row[-1] == True:
                continue

            if verbose:
                print('Looking for incoming paths to', name)

            paths = []

            for m in temp_nodes['UniProtID'].values:
                if m == n:
                    pass

                else:
                    paths.append(self.pathfinder(m, n))
                   
            lens = np.array([len(i) for i in paths])

            aux = []
            top = 1
            pending = True

            while pending:
                if top in lens:
                    if verbose:
                        print('* Found %d shortest paths to query' %np.sum(
                                                                   lens == top))
                        print('with %d intermediates' %top)

                    aux = np.unique([paths[i] for i in
                                     np.where(lens == top)[0]])
                    pending = False
                   
                else:
                    top += 1

            to_add.update(aux)

        if verbose:
            print('Adding %d intermediate nodes to the query' %len(to_add))

        self.update_query(to_add)

    def update_tables(self):
        '''
        Updates the nodes and edges information tables according to the current
        state of the network.

        This function creates/modifies the following attributes:
            - nodes, edges.
        '''

        # Nodes information:
        self.nodes = pd.DataFrame(index=range(self.network.vcount()),
                                  columns=['Name', 'UniProtID', 'Indegree',
                                           'Outdegree', 'Input'])
        self.nodes['Name'] = self.network.vs['label']
        self.nodes['UniProtID'] = self.network.vs['name']
        self.nodes['Indegree'] = self.network.vs.indegree()
        self.nodes['Outdegree'] = self.network.vs.outdegree()

        if self.mode == 'genesymbol':
            self.nodes['Input'] = [i in self.inputs for i in self.nodes['Name']]

        elif self.mode == 'uniprot':
            self.nodes['Input'] = [i in self.inputs for i in
                                   self.nodes['UniProtID']]

        else:
            raise Exception('ERROR: Invalid mode selection. ' +
                            'Please use "genesymbol" or "uniprot"')

        # Edges information:
        self.edges = pd.DataFrame(index=range(self.network.ecount()),
                                  columns=['Pair', 'Name', 'Refs',
                                           'Interaction (+, -)', 'Type'])

        # Source/Target node pairs as node IDs
        self.edges['Pair'] = [(e.source, e.target) for e in self.network.es]

        # Source/Target node pairs as node names
        self.edges['Name'] = [(self.nodes['Name'][i], self.nodes['Name'][j])
                              for (i, j) in  self.edges['Pair']]

        # Number of references supporting that edge
        self.edges['Refs'] = [len(self.network.es[i]['references']) for i in
                              self.edges.index]

        # Type of interaction (1: Activation | -1: Inhibition)
        aux = [(e['dirs'].is_stimulation(), e['dirs'].is_inhibition()) for e in
               self.network.es]
        self.edges['Interaction (+, -)'] = aux
        self.edges['Type'] = [('-' if i == (False, True) else '+') for i in aux]

    def reduce_inout(self, max_in=5, max_out=5, breakif=5, verbose=False):
        '''
        Takes a network and reduces the number of edges according to the
        maximum in/out degree given by the arguments. Edges are removed
        according to the increasing number of references (ones with less
        references are erased first). Also removes all incoming edges to any
        input node.

        This function modifies the following attributes:
            - nodes, edges, network.

        * Arguments:
            - max_in [int]: Optional, 5 by default. Determines the maximum
              number of incoming edges for a node in order to reduce the
              network complexity.
            - max_out [int]: Optional, 5 by default. Determines the maximum
              number of outgoing edges for a node in order to reduce network
              complexity.
            - breakif [int]: Optional, 5 by default. Sets a threshold for which
              edges whose number of references are equal or above, they cannot
              be removed.
            - verbose [bool]: Optional, False by default. Sets on/off the
              printing of messages.
        '''

        # Number incoming/outgoing edges for each node in the network
        ind = self.network.indegree()
        outd = self.network.outdegree()

        # List of node IDs whose in/out degree is over maximum
        high_ind = [i for i in range(self.network.vcount()) if ind[i] > max_in]
        high_outd = [i for i in range(self.network.vcount()) if outd[i] >
                     max_out]
       
        to_remove = set()
       
        # Iterate over high indegree nodes (n)
        for n in high_ind:
            if verbose:
                print('{} has {} incoming edges.'.format(self.nodes['Name'][n],
                                                         ind[n]))

            # Incoming edges to n
            in_eds = self.edges.loc[[i for i in self.edges.index if
                                    self.edges.loc[i, 'Pair'][1] == n],
                                    :].sort_values(by='Refs')
            i = 0
            # Remove incoming edges until threshold is reached
            while len(in_eds) > max_in:
                # Be sure that both nodes implied in the edge are not getting
                # disconnected from the network
                if (outd[in_eds.iloc[i]['Pair'][0]] > 1 and
                    ind[in_eds.iloc[i]['Pair'][1]] > 1):
                    # If the number of supporting references is above breakif
                    # threshold, stop removing edges
                    if in_eds.iloc[i]['Refs'] >= breakif:
                        if verbose:
                            print('* No further edges removed,',)
                            print('supporting references =>', breakif)
                            break

                    elif verbose:
                        print('  Removing edge {} --{}--> {}'.format(
                                                     in_eds.iloc[i]['Name'][0],
                                                     in_eds.iloc[i]['Refs'],
                                                     in_eds.iloc[i]['Name'][1]))

                    outd[in_eds.iloc[i]['Pair'][0]] -= 1
                    ind[in_eds.iloc[i]['Pair'][1]] -= 1
                    to_remove.update([in_eds.index[i]])

                    self.edges.drop(in_eds.index[i], inplace=True)
                    in_eds.drop(in_eds.index[i], inplace=True)

                else:
                    i += 1
                   
        # Iterate over high outdegree nodes (n)
        for n in high_outd:
            if verbose:
                print('{} has {} outcoming edges.'.format(
                                                 self.nodes['Name'][n], outd[n]))

            # Outgoing edges from n
            out_eds = self.edges.loc[[i for i in self.edges.index if
                                     self.edges.loc[i, 'Pair'][0] == n],
                                     :].sort_values(by='Refs')
            i = 0
            # Remove outgoing edges until threshold is reached
            while len(out_eds) > max_out:
                # Be sure that both nodes implied in the edge are not getting
                # disconnected from the network
                if (outd[out_eds.iloc[i]['Pair'][0]] > 1 and
                    ind[out_eds.iloc[i]['Pair'][1]] > 1):
                    # If the number of supporting references is above breakif
                    # threshold, stop removing edges
                    if out_eds.iloc[i]['Refs'] >= breakif:
                        if verbose:
                            print('* No further edges removed,',)
                            print('supporting references =>', breakif)
                            break

                    if verbose:
                        print('  Removing edge {} --{}--> {}'.format(
                                                    out_eds.iloc[i]['Name'][0],
                                                    out_eds.iloc[i]['Refs'],
                                                    out_eds.iloc[i]['Name'][1]))

                    outd[out_eds.iloc[i]['Pair'][0]] -= 1
                    ind[out_eds.iloc[i]['Pair'][1]] -= 1
                    to_remove.update([out_eds.index[i]])
                    self.edges.drop(out_eds.index[i], inplace=True)
                    out_eds.drop(out_eds.index[i], inplace=True)
                   
                else:
                    i += 1
        # Let us remove all incoming edges of the input nodes
        aux = [i[0] for i in self.edges.iterrows() if i[1][0][1] in
               self.nodes[self.nodes['Input']].index.values]
        to_remove.update(aux)

        # Update network and nodes information
        self.network.delete_edges(to_remove)
        self.update_tables()

    def pathfinder(self, a, b):
        '''
        Get list of intermediate nodes of the shortest path connecting a and b.

        * Arguments:
            - a [str]: Name of the origin node.
            - b [str]: Name of the destination node.

        * Returns:
            - [list]: Names [str] of the intermediate nodes of the shortest
              path from a to b (not included).
        '''

        path = self.pa.dgraph.get_shortest_paths(a, b)[0][1:-1]

        if self.mode == 'genesymbol':
            return [self.pa.dgraph.vs[i]['label'] for i in path]

        elif self.mode == 'uniprot':
            return [self.pa.dgraph.vs[i]['name'] for i in path]

        else:
            raise Exception('ERROR: Invalid mode selection. ' +
                            'Please use "genesymbol" or "uniprot"')

    def update_query(self, nodes):
        '''
        Updates the query dictionary containing the protein names and their
        corresponding UniProt IDs with the nodes passed in the argument.

        This function modifies the following attributes:
            - query.

        * Arguments:
            - nodes [list]: Contains all the node names [str] for which the
              UniProt IDs are to be retrieved. Apart from list, any other type
              of iterable can be used as long as its elements are [str] (e.g.:
              a set, a tuple, ...).
        '''

        if self.mode == 'genesymbol':
            self.query.update([(str(n), self.mappr.map_name(n, 'genesymbol',
                                                            'uniprot')[0])
                               for n in nodes if str(n) not in
                               self.query.values()])

        elif self.mode == 'uniprot':
            self.query.update([(self.mappr.map_name(n, 'uniprot',
                                                    'genesymbol')[0], str(n))
                               for n in nodes if str(n) not in
                               self.query.values()])

        else:
            raise Exception('ERROR: Invalid mode selection. ' +
                            'Please use "genesymbol" or "uniprot"')

    def to_sif(self, filename):
        '''
        Writes a .sif file containing all the edges from the network.

        * Arguments:
            - filename [str]: File name or path under which the file is to be
              saved. The .sif file extension can be obviated (in such case, it
              will be added automatically).
        '''
        if filename[-4:] != '.sif':
            filename += '.sif'

        interact = [(-1 if i == '-' else 1) for i in self.edges.loc[:, 'Type']]
        lines = ['{}\t{}\t{}'.format(self.edges.loc[i, 'Name'][0],
                                     interact[i],
                                     self.edges.loc[i, 'Name'][1])
                 for i in self.edges.index.values]

        with open(filename, 'w') as f:
            f.write('\n'.join(lines))

    def plot_adjacency(self, filename=None):
        '''
        Plots the adjacency matrix of the network. If a file name or path is
        passed, also saves the figure.

        * Arguments:
            - filename [str]: Optional, None by default. File name or path
              where the figure is to be stored.
        '''

        fig, ax = plt.subplots(figsize=(10, 10))

        adj = np.array(self.network.get_adjacency().data)
        names = self.network.vs()['label']
        rng = range(len(names))

        ax.imshow(adj, cmap='gray')
        ax.set_xlabel('Target')
        ax.set_ylabel('Source')
        ax.set_xticks(rng)
        ax.set_xticklabels(names, rotation=90, fontsize=9)
        ax.set_yticks(rng)
        ax.set_yticklabels(names, fontsize=9)

        fig.tight_layout()

        if filename:
            fig.savefig(filename)

################################################################################

class PathFinder(object): # TODO: Add/re-implement network reduction?
    '''
    Main class. Contains a copy of PyPath and ultimately the PKN built from the
    nodes set by the user. The latter is a directed graph object from igraph.
    Futher information about PyPath or igraph can be found in their respective
    documentation.

    * Attributes:
        - pa [pypath.main.PyPath]: An initialized PyPath instance (see the
          package reference for further details).
        - mappr [pypath.mapping.Mapper]: PyPath mapping object. Helps
          converting the provided identifiers to UniProt IDs according to mode.
        - inputs [set]: Contains the UniProt IDs [str] of the nodes that are
          to be considered as inputs. These will not have incoming edges in the
          final final PKN.
        - nodes [set]: Contains the UniProt IDs [str] of all the nodes which
          are not considered as inputs. These will have at least one incoming
          edge in the final PKN.
        - node_table [pandas.DataFrame]: Contains information about the nodes
          in the network. More specifically: 'Name', 'UniProtID', 'Indegree',
          'Outdegree' and 'Input' (True if it is, False otherwise).
        - edge_table [pandas.DataFrame]: Contains some information about the
          edges in the network, like: 'Pair' (source-target nodes as indices),
          'Name' (of the node pair), 'Refs' (number of), 'Interaction (+, -)'
          (indicating if there is supporting reference for the edge being
          activatory (+) or inhibitory (-)) and 'Type' (considered inhibitory
          if and only if the reference supports that relation). Please note
          that the last information should be confirmed manually in the
          literature as it can be ambiguous, incomplete or incorrect (See
          'Interaction (+, -)' column).
        - pkn [igraph.Graph]: Contains the PKN as an igraph object. See the
          corresponding documentation for more information about object methods
          and attributes.
        - mode [str]: Specifies the type of identifiers used on the input.
          Available options are: 'uniprot', 'genesymbol', 'entrez', 'refseq',
          'ensp', 'enst', 'ensg', 'hgnc', 'gi', 'embl' and 'embl_id'.
    '''
    def __init__(self, pa=None, filename=None):
        '''
        Initializes the pkn graph and loads the ID mapper. Creates the PyPath
        instance and downloads all the data from scratch if not loaded from a
        file or another PyPath instance. If using a PyPath object already
        loaded in the current environment, the following does not apply.
        When loading PyPath, edges with only high-throughput references are
        removed. Also loads the gene symbols mapping to UniProt IDs, pathway
        attributes and edge directions. If the instance is loaded from scratch
        (filename and pa not provided), user will be prompted to input a new
        file name/path to save the new pickle (if none is provided, it does
        nothing).

        This function creates the following attributes:
            - pkn, pa, mappr.

        * Arguments:
            - pa [pypath.main.PyPath]: Optional, None by default. An
              initialized PyPath instance (see the package reference for
              further details). If none is provided, checks if filename is
              provided.
            - filename [str]: Optional, None by default. File name or path to
              the file where the PyPath pickle should be loaded. If none is
              provided then all the data is downloaded from scratch.
        '''

        self.pkn = igraph.Graph(directed=True)
        self.mappr = pypath.mapping.Mapper(9606)

        if pa:
            self.pa = pa

        else:
            self.pa = pypath.PyPath()

            if filename:
                self.pa.init_network(pfile=filename)

            else:
                self.pa.init_network()

            self.pa.remove_htp()
            self.pa.genesymbol_labels()
            self.pa.pathway_attributes()
            self.pa.load_all_pathways()
            self.pa.get_directed()

            if not filename:
                filename = raw_input('Enter name/path to save the pickle: ')
                if filename is not '':
                    self.pa.save_network(filename)
                    print('PyPath instance saved at %s' %filename)

    def init_nodes(self, inputs, nodes, mode='genesymbol', check=True,
                   verbose=False):
        '''
        Loads the nodes provided by the user into the PKN. Nodes should be
        provided in two different lists, Namely inputs and nodes. Several
        types of identifiers are supported and should be specified in the
        keyword argument mode.

        This function creates/modifies the following attributes:
            - mode, inputs, nodes, pkn.

        * Arguments:
            - inputs [list]: Contains the IDs [str] of the nodes that are to be
              considered as inputs. These will not have incoming edges in the
              final final PKN.
            - nodes [list]: Contains the IDs [str] of all the nodes which are
              not considered as inputs. These will have at least one incoming
              edge in the final PKN.
            - mode [str]: Optional, 'genesymbol' by default. Specifies the type
              of identifiers used on the input. Other available options are:
              'uniprot', 'entrez', 'refseq', 'ensp', 'enst', 'ensg', 'hgnc',
              'gi', 'embl' and 'embl_id'.
            - check [bool]: Optional, True by default. Wether to check the IDs
              if already provided as 'uniprot' identifiers (can save some time,
              specially for large queries).
            - verbose [bool]: Optional, False by default. Sets on/off the
              printing of messages.
        '''

        self.mode = mode

        if mode == 'uniprot' and not check:
            self.inputs = set(inputs)
            self.nodes = set(nodes)

        else: # Converting IDs to UniProt / Checking UniProt IDs
            if verbose:
                msg = ('Checking UniProt IDs' if self.mode is 'uniprot'
                       else 'Converting %s to UniProt IDs' %self.mode)
                print(msg)
            self.inputs = set(self.convert(inputs))
            self.nodes = set(self.convert(nodes))

        # Adding inputs to the pkn
        for v in self.inputs:
            self.add_node(v, is_input=True)

        # Adding other nodes to the pkn
        for v in self.nodes:
            if v not in self.pkn.vs['name']:
                self.add_node(v)

    def build_network(self, verbose=False):
        '''
        For each node in the query, checks if there is already any direct
        target in the set. In such case, that edge is added to the PKN. Then
        the connectivity of each node is checked, so that input nodes have at
        least an out-degree over 0 (keeping in-degree always 0) and the other
        nodes an in-degree over 0. If a node does not fulfil the condition, the
        sorthest paths to any other node in the query are cheked. Among these,
        the shortest are added to the PKN. In the case that no path is to be
        found for a specific node, a warning message is prited.

        This function modifies the following attributes:
            - pkn, node_table, edge_tabe.

        * Arguments:
            - verbose [bool]: Optional, False by default. Sets on/off the
              printing of messages.
        '''

        # First we build the direct connections
        for n in self.pkn.vs:
            out_conect, in_connect = [], []

            out_connect = [self.pa.dgraph.vs[i]['name'] for i in
                           self.pa.dgraph.neighbors(n['name'], mode='OUT')]
            out_connect = [i for i in out_connect if (i in self.pkn.vs['name']
                           and not self.pkn.vs.find(i)['is_input'])]

            if not n['is_input']:
                in_connect = [self.pa.dgraph.vs[i]['name'] for i in
                              self.pa.dgraph.neighbors(n['name'], mode='IN')]
                in_connect = [i for i in in_connect if i in self.pkn.vs['name']]

            ncon = len(in_connect) + len(out_connect)

            if ncon:
                if verbose:
                    print('%s has %d connection(s) to query' %(n['label'], ncon))
                if out_connect:
                    for m in out_connect:
                        self.add_edge(n, self.pkn.vs.find(m))

                if in_connect:
                    for m in in_connect:
                        self.add_edge(self.pkn.vs.find(m), n)


        # Now we connect the other nodes (if possible)
        for n in self.pkn.vs:
            if n['is_input']: # Input nodes must have indegree = 0 < outdegree
                connectivity = n.outdegree()

            else: # Non-input nodes must have indegree > 0
                connectivity = n.indegree()

            if connectivity < 1:
                if verbose:
                    print('* %s has no direct connections in query' % n['label'])

                # List of lists containing the vertices of shortest paths to
                # all query nodes
                paths = []

                for m in self.pkn.vs:
                    if m == n:
                        pass

                    else:
                        if n['is_input'] and m['is_input']:
                            pass

                        elif n['is_input'] and not m['is_input']:
                            paths.extend(self.pathfinder(n, m))

                        else:
                            paths.extend(self.pathfinder(m, n))

                lens = np.array([len(p) for p in paths])

                shortest = []
                top = 2
                pending = True

                if sum(lens) == 0:
                    pending = False
                    print('*** WARNING ***: No connection found for', n['label'])

                while pending:
                    if top in lens:
                        aux = np.where(lens == top)[0]
                        if verbose:
                            print('\tFound %d shortest paths to' %len(aux),)
                            print('query with %d intermediates' %(top - 2))

                        shortest = [paths[i] for i in aux]
                        pending = False

                    else:
                        top += 1

                for path in shortest:
                    if verbose:
                        print('\t' + ' -> '.join([p['label'] for p in path]))

                    pairs = zip(path, path[1:])
                    for p0, p1 in pairs:
                        if p1['name'] not in self.pkn.vs['name']:
                            if verbose:
                                print('\t* Added node %s' %(p1['label']))

                            self.add_node(p1['name'])

                        self.add_edge(self.pkn.vs.find(p0['name']),
                                      self.pkn.vs.find(p1['name']))

        self.update_tables()

    def add_node(self, n, is_input=False):
        '''
        Adds a node to the PKN. All attributes of that node in PyPath are
        copied alongside.

        This function modifies the following attributes:
            - pkn.

        * Arguments:
            - n [str]: The UniProt ID of the node to be added.
            - is_input [bool]: Optional, False by default. Specifies if the
              newly created node has to be considered as an input.
        '''

        self.pkn.add_vertex(n)
        self.pkn.vs.find(n).update_attributes(self.pa.duniprot(n).attributes())
        self.pkn.vs.find(n).update_attributes({'is_input':is_input})

    def add_edge(self, a, b):
        '''
        Adds an edge to the PKN except if that edge already exists or the
        target node is an input. All attributes of that edge in PyPath are
        copied alongside.

        This function modifies the following attributes:
            - pkn.

        * Attributes:
            - a [igraph.Vertex]: Source vertex of the edge to be added. NOTE:
              this must be drawn from the PKN and not from PyPath.
            - b [igraph.Vertex]: Target vertex of the edge to be added. NOTE:
              this must be drawn from the PKN and not from PyPath.
        '''

        if b['is_input'] or bool(self.pkn.es.select(_source=a.index,
                                                    _target=b.index)):
            pass # Edge already exists in PKN or target is an input

        else:
            self.pkn.add_edge(a.index, b.index)
            self_eid = self.pkn.get_eid(a.index, b.index)
            pa_eid = self.pa.dgraph.get_eid(self.pa.duniprot(a['name']).index,
                                            self.pa.duniprot(b['name']).index)
            attrs = self.pa.dgraph.es[pa_eid].attributes()
            self.pkn.es[self_eid].update_attributes(attrs)

    def pathfinder(self, a, b): # XXX: EDIT DOCS
        '''
        Retrieves the shortest path between two given nodes.

        * Arguments:
            - a [igraph.Vertex]: Source node for the path search.
            - b [igraph.Vertex]: Target node for the path search.

        * Returns:
            - [list]: List of lists (paths) containing the nodes
              [igraph.Vertex] inolved on the shortest paths from a
              to b (both included). NOTE: nodes are drawn from PyPath instance.
        '''

        paths = self.pa.dgraph.get_shortest_paths(a['name'], b['name'])

        # Check that none of the paths has an input node in the middle
        npaths = [[self.pa.dgraph.vs[i] for i in path] for path in paths]

        for npath in npaths:
            is_input = [n['name'] in self.inputs for n in npath]
            if True in is_input:
                npaths.remove(npath) # In such case, remove it

        return npaths

    def convert(self, ids):
        '''
        Converts a given list of identifiers (whose type is specified in the
        attribute mode) into a list of UniProt IDs. If none of them can be
        mapped, an exception is raised. In case that only some of them could
        not be mapped, a warning message is printed and only the mapped IDs are
        returned.

        * Arguments:
            - ids [list]: Contains the identifiers [str] to be mapped.

        * Returns:
            - [list]: Mapped UniProt IDs [str].
        '''

        try: # Silent mode
           
            with silent():

                result = self.mappr.map_names(ids, self.mode, 'uniprot')

            if len(result) == 0:
                raise Exception('No IDs could be mapped. Invalid mode or IDs')

            elif len(result) < len(ids):
                print('*** WARNING ***: Some of the IDs could not be mapped')
                # TODO: Could also print(which IDs couldn't be mapped.)

        except Exception as err:
            raise err

        else:
            return result

    def update_tables(self):
        '''
        Updates the nodes and edges information tables according to the current
        state of the PKN.

        This function creates/modifies the following attributes:
            - node_table, edge_table.
        '''

        # Nodes information:
        self.node_table = pd.DataFrame(index=range(self.pkn.vcount()),
                                       columns=['Name', 'UniProtID',
                                                'Indegree', 'Outdegree',
                                                'Input'])
        self.node_table['Name'] = self.pkn.vs['label']
        self.node_table['UniProtID'] = self.pkn.vs['name']
        self.node_table['Indegree'] = self.pkn.vs.indegree()
        self.node_table['Outdegree'] = self.pkn.vs.outdegree()
        self.node_table['Input'] = self.pkn.vs['is_input']

        # Edges information:
        self.edge_table = pd.DataFrame(index=range(self.pkn.ecount()),
                                       columns=['Pair', 'Name', 'Refs',
                                                'Interaction (+, -)', 'Type'])

        # Source/Target node pairs as node IDs
        self.edge_table['Pair'] = [(e.source, e.target) for e in self.pkn.es]

        # Source/Target node pairs as node names
        self.edge_table['Name'] = [(self.pkn.vs[i]['label'],
                                    self.pkn.vs[j]['label'])
                                   for (i, j) in self.edge_table['Pair']]

        # Number of references supporting that edge
        self.edge_table['Refs'] = [len(e['references']) for e in self.pkn.es]

        # Type of interaction (1: Activation | -1: Inhibition)
        aux = [(e['dirs'].is_stimulation(), e['dirs'].is_inhibition())
               for e in self.pkn.es]

        self.edge_table['Interaction (+, -)'] = aux
        self.edge_table['Type'] = [('-' if i == (False, True) else '+')
                        for i in aux]

    def plot_adjacency(self, filename=None, figsize=None):
        '''
        Plots the adjacency matrix of the network. If a file name or path is
        passed, also saves the figure.

        * Arguments:
            - filename [str]: Optional, None by default. File name or path
              where the figure is to be stored.
            - figsize [tuple]: Optional, None by default (matplotlib default
              size). Specifies the figure size in inches (width, height).
        '''

        fig, ax = plt.subplots(figsize=figsize)

        adj = np.array(self.pkn.get_adjacency().data)
        names = self.pkn.vs()['label']
        rng = range(len(names))

        ax.imshow(adj, cmap='gray')
        ax.set_xlabel('Target')
        ax.set_ylabel('Source')
        ax.set_xticks(rng)
        ax.set_xticklabels(names, rotation=90, fontsize=9)
        ax.set_yticks(rng)
        ax.set_yticklabels(names, fontsize=9)

        fig.tight_layout()

        if filename:
            fig.savefig(filename)

    def to_sif(self, filename):
        '''
        Writes a .sif file containing all the edges from the network.

        * Arguments:
            - filename [str]: File name or path under which the file is to be
              saved. The .sif file extension can be obviated (in such case, it
              will be added automatically).
        '''

        if filename[-4:] != '.sif':
            filename += '.sif'

        interact = [(-1 if i == '-' else 1)
                    for i in self.edge_table.loc[:, 'Type']]
        lines = ['{}\t{}\t{}'.format(self.edge_table.loc[i, 'Name'][0],
                                     interact[i],
                                     self.edge_table.loc[i, 'Name'][1])
                 for i in self.edge_table.index.values]

        with open(filename, 'w') as f:
            f.write('\n'.join(lines))

################################################################################

#--------------------------------- FUNCTIONS ----------------------------------#
def convert_names(filename, keys):
    '''
    Function that reads a given .sif file and substitutes all names specified
    on the passed dictionary.

    * Arguments:
        - filename [str]: File name or path under for the source .sif file.
        - keys [dict]: Key/value pairs [str] containing the new and old names
          respectively. For instance, if we only want to rename a node named
          'A' as 'a' then keys = {'a':'A'}.
    '''

    def subs(string):
        for k, v in iteritems(keys):
            string = string.replace(v, k)

        return string

    with open(filename) as f:
        lines = [r for r in f]

    with open(filename, 'w') as f:
        f.write(''.join(map(subs, lines)))

def read_sif(filename):
    '''
    Reads a .sif file and returns a list containing all the nodes in it.

    * Arguments:
        - filename [str]: File name or path under for the source .sif file.

    * Returns:
        - [list]: Contains all the unique node names [str] contained on the
          file.
    '''

    nodes = set()
    with open(filename) as f:
        for l in f.readlines():
            # Ignore comments, remove newline characters and split by separator
            if '\t' in l:
                aux = l.split('#')[0].strip().split('\t')

            else:
                aux = l.split('#')[0].strip().split(' ')

            if len(aux) is 3:
                nodes.update([aux[0], aux[2]])

    return list(nodes)

#------------------------------------------------------------------------------#
