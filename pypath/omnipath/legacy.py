#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright 2014-2023
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: see the file `README.rst`
#  Contact: Dénes Türei (turei.denes@gmail.com)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      https://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: https://pypath.omnipathdb.org/
#

import os
import copy

from pypath.legacy import main
from pypath.resources import data_formats
from pypath.core import annot
from pypath.core import intercell
from pypath.core import complex
from pypath.core import ptm
from pypath.share import settings
from pypath.share import session as session_mod


class OmniPath(session_mod.Logger):
    
    
    def __init__(
        self,
        output_dir = None,
        network_pickle = None,
        annotation_pickle = None,
        intercell_pickle = None,
        complex_pickle = None,
        enz_sub_pickle = None,
        load_network = True,
        load_complexes = True,
        load_annotations = True,
        load_intercell = True,
        load_enz_sub = True,
    ):
        
        if not hasattr(self, '_log_name'):
            
            session_mod.Logger.__init__(self, name = 'omnipath')
        
        self.output_dir = output_dir
        self.network_pickle = network_pickle
        self.annotation_pickle = annotation_pickle
        self.intercell_pickle = intercell_pickle
        self.complex_pickle = complex_pickle
        self.enz_sub_pickle = enz_sub_pickle
        
        self.do_load_network = load_network
        self.do_load_complexes = (
            load_complexes or
            load_annotations or
            load_intercell
        )
        self.do_load_annotations = load_annotations or load_intercell
        self.do_load_intercell = load_intercell
        self.do_load_enz_sub = load_enz_sub
        
        self.main()
    
    
    def main(self):
        
        self.load()
    
    
    def load(self):
        
        self.load_complex()
        self.load_network()
        self.load_annotations()
        self.load_intercell()
        self.load_enz_sub()
    
    
    def load_complex(self):
        
        if not self.do_load_complexes:
            
            return
        
        self.complex = complex.get_db(
            pickle_file = self.ensure_path_exists(self.complex_pickle)
        )
    
    
    def load_network(self):
        
        if not self.do_load_network:
            
            return
        
        self.network = main.PyPath()
        
        if os.path.exists(self.network_pickle):
            
            self.network.init_network(pfile = self.network_pickle)
            
        else:
            
            network_input = copy.deepcopy(data_formats.omnipath)
            network_input.update(data_formats.ligand_receptor)
            network_input.update(data_formats.ptm_misc)
            self.network.load_omnipath(
                omnipath = network_input,
                remove_htp = False,
            )
    
    
    def load_annotations(self):
        
        if not self.do_load_annotations:
            
            return
        
        self.annot = annot.get_db(
            pickle_file = self.ensure_path_exists(self.annotation_pickle)
        )
    
    
    def load_intercell(self):
        
        if not self.do_load_intercell:
            
            return
        
        self.intercell = (
            intercell.get_db(
                pickle_file = self.ensure_path_exists(self.intercell_pickle)
            )
        )
    
    
    def load_enz_sub(self):
        
        if not self.do_load_enz_sub:
            
            return
        
        self.enz_sub = (
            ptm.get_db(
                pickle_file = self.ensure_path_exists(self.enz_sub_pickle)
            )
        )
    
    
    @staticmethod
    def ensure_path_exists(path):
        
        return path if path and os.path.exists(path) else None
