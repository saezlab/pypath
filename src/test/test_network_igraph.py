#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2018
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
Tests for the `igraph` based implementation of pypath's network building
tools. Currently in `pypath.main.PyPath`.
"""

import pytest

import pypath.main as main
import pypath.settings as settings
import pypath.data_formats as data_formats


class TestPyPath(object):
    
    
    def test_complex_expansion(self):
        
        input_param = {'Signor': data_formats.pathway['signor']}
        
        # a network with all complexes expanded to single proteins
        settings.setup(network_complex_expansion = True)
        netw_complex_expanded = main.PyPath()
        netw_complex_expanded.init_network(input_param)
        
        # a network with complexes as entities (vertices) in it
        settings.setup(network_complex_expansion = False)
        netw_complex_entities = main.PyPath()
        netw_complex_entities.init_network(input_param)
        
        # collecting all interactions of the complexes
        expanded_connections = set()
        
        for v in netw_complex_entities.graph.vs:
            
            if v['type'] != 'complex':
                
                continue
            
            for vid_nb in netw_complex_entities.graph.neighbors(v.index):
                
                for cplex_comp in v['name'].components.keys():
                    
                    expanded_connections.add((
                        netw_complex_entities.vs['name'][vid_nb],
                        cplex_comp,
                    ))
        
        # checking if all these connections exist in the expanded network
        missing = set()
        
        for uniprot1, uniprot2 in expanded_connections:
            
            if not netw_complex_expanded.get_edge(uniprot1, uniprot2):
                
                missing.add(uniprot1, uniprot2)
        
        assert not missing
