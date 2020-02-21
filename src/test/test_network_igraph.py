#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2020
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#                  Nicolàs Palacio
#                  Olga Ivanova
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

import pypath.resources.data_formats as data_formats
from pypath.legacy import main
import pypath.share.settings as settings


class TestPyPath(object):

    def test_complex_expansion(self):
        input_param = {'Signor': data_formats.pathway['signor']}

        settings.setup(network_pickle_cache=False)
        # a network with all complexes expanded to single proteins
        settings.setup(network_expand_complexes=True)
        netw_complex_expanded = main.PyPath()
        netw_complex_expanded.keep_raw = True
        netw_complex_expanded.init_network(input_param)

        # a network with complexes as entities (vertices) in it
        settings.setup(network_expand_complexes=False)
        netw_complex_entities = main.PyPath()
        netw_complex_entities.init_network(input_param)

        # collecting all interactions of the complexes
        expanded_connections = set()

        for v in netw_complex_entities.graph.vs:
            if v['type'] != 'complex':
                continue

            for vid_nb in netw_complex_entities.graph.neighbors(v.index):
                names_nb = netw_complex_entities.graph.vs[vid_nb]['name']
                # if the neighbor is another complex
                uniprots_nb = (
                    names_nb.components.keys()
                    if hasattr(names_nb, 'components') else
                    (names_nb,)
                )

                for uniprot_nb in uniprots_nb:
                    for uniprot_comp in v['name'].components.keys():
                        expanded_connections.add(tuple(sorted((
                            uniprot_nb,
                            uniprot_comp,
                        ))))

        # checking if all these connections exist in the expanded network
        missing = set()

        for uniprot1, uniprot2 in expanded_connections:
            if (
                    not netw_complex_expanded.get_edge(uniprot1, uniprot2) and
                    uniprot1 != uniprot2
            ):
                missing.add((uniprot1, uniprot2))

        assert not missing

    def test_orthology_translation(self):
        """
        We test this on a nice complicated example case when 3 human
        genes map to 4 mouse orthologues. We select one interacting partner
        and check if it has one and only one interaction with this partner
        with the correct direction and sign.
        """

        pa = main.PyPath()
        pa.init_network({'s': data_formats.pathway['signor']})
        pa.homology_translation(target = 10090)
        
        adcy1_vs = pa.graph.vs.select(name = 'O88444')
        gnas_1_vs = pa.graph.vs.select(name = 'P63094')
        gnas_2_vs = pa.graph.vs.select(name = 'Q6R0H7')
        gnas_3_vs = pa.graph.vs.select(name = 'Q6R0H6')
        gnas_4_vs = pa.graph.vs.select(name = 'Q9Z0F1')
        
        assert all(
            len(vs) == 1
            for vs in (adcy1_vs, gnas_1_vs, gnas_2_vs, gnas_3_vs, gnas_4_vs)
        )
        
        adcy1_i = adcy1_vs[0].index
        
        for gnas_vs in (gnas_1_vs, gnas_2_vs, gnas_3_vs, gnas_4_vs):
            
            gnas_i = gnas_vs[0].index
            
            es = pa.graph.es.select(_between = ((gnas_i,), (adcy1_i,)))
            
            assert len(es) == 1
            
            edge = es[0]
            
            assert edge['dirs'].which_dirs()[0][1] == 'O88444'
            assert (
                edge['dirs'].get_sign(
                    *edge['dirs'].which_dirs()
                ) == [True, False]
            )
