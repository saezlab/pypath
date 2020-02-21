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
#                  Olga Ivanova
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

import pypath.resources.network as netres
from pypath import main


class TestInteraction(object):

    def test_contains(self):
        sr = netres.pathway['signor']
        sk = sr.key
        si = 'SIGNOR'

        pa = main.PyPath()
        pa.init_network(netres.pathway)

        assert (
                [sr in a for a in pa.graph.es['attrs']] ==
                [sk in a for a in pa.graph.es['attrs']] ==
                [si in a for a in pa.graph.es['attrs']]
        )

    def test_filter_via(self):
        pa = main.PyPath()
        pa.init_network(netres.ptm_misc)

        assert any(
            len(list(a.evidences.filter(via='MIMP')))
            for a in pa.graph.es['attrs']
        )
