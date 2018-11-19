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


import sys
import pytest

import pypath.dataio as dataio
import pypath.data_formats as data_formats
import pypath.settings as settings


exclude = {
    'reactome_interactions'
}
network_methods = []

for var in data_formats.__dir__():
    
    var = getattr(data_formats, var)
    
    if isinstance(var, dict):
        
        for k, v in var.items():
            
            if hasattr(v, 'inFile'):
                
                method_name = getattr(v, 'inFile')
                
                if (
                    hasattr(dataio, method_name) and
                    method_name not in exclude
                ):
                    
                    network_methods.append((method_name, v.inputArgs))


class TestDataio(object):
    
    @pytest.mark.parametrize('method, args', network_methods)
    def test_network_data_methods(self, method, args, cachedir):
        
        settings.setup(cachedir = cachedir)
        
        method = getattr(dataio, method)
        result = method(**args)
        
        assert len(result)
