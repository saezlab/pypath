#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#
#  Copyright (c) 2014-2015 - EMBL-EBI
#
#  File author(s): Dénes Türei (denes@ebi.ac.uk)
#
#  Distributed under the GNU GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://www.ebi.ac.uk/~denes
#

import bioigraph

net = bioigraph.BioGraph()
net.init_network(exclude = ['intact', 'acsn', 'reactome', 'nci-pid'])
net.remove_htp()
net.third_source_directions()
net.save_network('cache/default_network.pickle')