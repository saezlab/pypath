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

from future.utils import iteritems

import imp
import collections

import pypath.session_mod as session_mod
import pypath.main as main_mod
import pypath.network as network_mod
import pypath.data_formats as data_formats
import pypath.annot as annot
import pypath.intercell as intercell


CellPhoneDBGene = collections.namedtuple(
    'CellPhoneDBGene',
    [
        'gene_name',
        'uniprot',
        'hgnc_symbol',
        'ensembl',
    ],
)


CellPhoneDBProtein = collections.namedtuple(
    'CellPhoneDBProtein',
    [
        'uniprot',
        'adhesion',
        'cytoplasm',
        'entry_name',
        'extracellular',
        'integrin_interaction',
        'other',
        'other_desc',
        'pdb_id',
        'pdb_structure',
        'peripheral',
        'receptor',
        'receptor_desc',
        'secreted_desc',
        'secreted_highlight',
        'secretion',
        'stoichiometry',
        'tags_description',
        'transmembrane',
        'transporter',
        'tags',
        'tags_reason',
    ],
)


CellPhoneDBInteraction = collections.namedtuple(
    'CellPhoneDBInteraction',
    [
        'comments_interaction',
        'dlrp',
        'family',
        'iuphar',
        'multidata_name_1',
        'multidata_name_2',
        'score_1',
        'score_2',
        'source',
    ],
)


CellPhoneDBComplex = collections.namedtuple(
    'CellPhoneDBComplex',
    [
        'name',
        'uniprot_1',
        'uniprot_2',
        'uniprot_3',
        'uniprot_4',
        'adhesion',
        'cytoplasm',
        'extracellular',
        'integrin_interaction',
        'other',
        'other_desc',
        'peripheral',
        'receptor',
        'receptor_desc',
        'secreted_desc',
        'secreted_highlight',
        'secretion',
        'transmembrane',
        'transporter',
        'pdb_structure',
        'pdb_id',
        'stoichiometry',
        'comments_complex',
    ],
)


class CellphoneDB(session_mod.Logger):
    
    def __init__(
            self,
            network = None,
            annotation = None,
            network_param = None,
            annot_param = None,
            omnipath = False,
        ):
        
        session_mod.Logger.__init__(self, name = 'cellphonedb')
        
        self.network = network
        self.annotation = annotation
        self.network_param = network_param or {}
        self.annot_param = annot_param or {}
        self.omnipath = omnipath
    
    
    def reload(self):
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        import imp
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    
    def main(self):
        
        self.setup_network()
        self.setup_annotation()
        self.setup_complex()
    
    
    def setup_network(self):
        
        if self.network is None:
            
            if not self.omnipath and 'lst' not in self.network_param:
                
                self.network_param['lst'] = data_formats.pathway
            
            self.network = main_mod.PyPath()
            self.network.init_network(**self.network_param)
        
        if isinstance(self.network, main_mod.PyPath):
            
            self.network = network_mod.Network(records = self.network)
            
        if isinstance(self.network, network_mod.Network):
            
            pass
    
    
    def setup_annotation(self):
        
        if self.annotation is None:
            
            self.annotation = intercell.IntercellAnnotation(
                **self.annot_param
            )
