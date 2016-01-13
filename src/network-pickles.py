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

import pypath
from pypath.common import *
from pypath.data_formats import *

net = pypath.Pypath(ncbi_tax_id = 9606)

net.init_network()
net.save_network(pfile = 'cache/default_network.pickle')
net.load_resources(lst={'acsn': ugly['acsn']})
net.save_network(pfile = 'cache/default_plus_acsn.pickle')

net = pypath.Pypath(ncbi_tax_id = 9606)
net.init_network(exclude = ['intact'])
net.save_network(pfile = 'cache/default_network_wo-intact.pickle')

net.load_resources(lst={'acsn': ugly['acsn']})
net.save_network(pfile = 'cache/default_plus_acsn_wo-intact.pickle')

net = pypath.Pypath(ncbi_tax_id = 9606)
net.init_network(exclude = ['intact'])
net.remove_htp()
net.save_network(pfile = 'cache/default_network_wo-intact_ltp-only.pickle')

net.load_resources(lst={'acsn': ugly['acsn']})
net.save_network(pfile = 'cache/default_plus_acsn.pickle')

net.load_resources(lst={'li2012': ugly['li2012']})
net.load_resources(lst={'mimp': good['mimp']})
net.load_resources(lst={'pnetworks': good['pnetworks']})
net.load_resources(lst={'psite_noref': good['psite_noref']})
net.save_network(pfile = 'cache/default_plus_acsn_phospho.pickle')
