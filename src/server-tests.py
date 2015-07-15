#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#
#  Copyright (c) 2014-2015 - EMBL-EBI
#
#  File author(s): Dénes Türei (denes@ebi.ac.uk)
#
#  This code is not for public use.
#  Please do not redistribute.
#  For permission please contact me.
#
#  Website: http://www.ebi.ac.uk/~denes
#

import bioigraph
from bioigraph import chembl
from bioigraph.common import *
from bioigraph.data_formats import *
from bioigraph import server

mysql_gelati = (None,'mapping_gelati')
mysql_chembl = (None,'chembl_ebi')

net = bioigraph.BioGraph(9606, mysql=mysql_gelati, name="demo")

net.init_network(pfile = 'cache/plus_phospho.pickle')
#net.load_resources(lst={'mimp': good['mimp']})
#net.load_resources(lst={'pnetworks': good['pnetworks']})
#net.load_resources(lst={'psite_noref': good['psite_noref']})
#net.load_resources(lst={'acsn': ugly['acsn']})
#net.save_network(pfile = 'cache/plus_phospho.pickle')

server.Rest(net, 33333)
