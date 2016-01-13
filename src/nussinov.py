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

import pypath
import copy
from itertools import chain
from pypath.data_formats import best, good, ugly, transcription
from pypath import dataio

mysql_gelati = (None,'mapping_gelati')
mysql_chembl = (None,'chembl_ebi')

net = pypath.Pypath(9606)

net.init_network(pfile = 'cache/plus_phospho.pickle')
net.load_resources(transcription)
net.load_go()
net.go_dict()

# number of TF-target links:
len([e for e in net.graph.es if 'TF' in e['type']])
# number of PPIs:
len([e for e in net.graph.es if 'PPI' in e['type']])
# number of TFs:
len(net.transcription_factors())
# number of parallel TF & PPI links:
len([e for e in net.graph.es if len(e['sources_by_type']) == 2])
