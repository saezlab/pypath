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

import imp
import collections

import pypath.session_mod as session_mod
import pypath.annot as mod_annot


IntercellRole = collections.namedtuple(
    'IntercellRole',
    ['source', 'role'],
)


class IntercellAnnotation(session_mod.Logger):
    
    
    def __init__(
            self,
            annot = None,
            annot_args = None,
        ):
        
        session_mod.Logger.__init__(self, name = 'intercell')
        
        self.annot = annot
        self.annot_args = annot_args or {}
        self.set_annot()
    
    
    def reload(self):
        """
        Reloads the object from the module level.
        """
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    
    def set_annot(self):
        """
        Creates an ``annot.AnnotationTable`` object if not provided.
        """
        
        self.annot = (
            self.annot or
            mod_annot.AnnotationTable(
                keep_annotators = True,
                **self.annot_args,
            )
        )
    
    
    def collect_receptors(self):
        """
        Creates a consensus annotation of plasma membrane receptors.
        """
        
        self.receptors_by_resource = {}
        
        self.add_receptors_cellphonedb()
        self.add_receptors_hpmr()
        self.add_receptors_go()
        self.add_receptors_surfaceome()
        
        self.receptors = set.union(*self.receptors_by_resource.values())
    
    
    def collect_ecm(self):
        
        self.ecm_by_resource = {}
        
        self.add_ecm_go()
        self.add_ecm_matrisome()
        self.add_ecm_matrixdb()
        
        self.ecm_proteins = set.union(*self.ecm_by_resource.values())
    
    
    def add_receptors_hpmr(self):
        
        self.receptors_by_resource['HPMR'] = (
            self.annot.annots['HPMR'].to_set()
        )
    
    
    def add_receptors_cellphonedb(self):
        
        self.receptors_by_resource['CellPhoneDB'] = (
            self.annot.annots['CellPhoneDB'].get_subset(
                receptor = bool,
                transmembrane = True,
            )
        )
    
    
    def add_receptors_go(self):
        
        self.receptors_by_resource['GO'] = (
            self.annot.annots['GO_Intercell'].get_subset(
                mainclass = 'receptors',
            )
        )
    
    
    def add_receptors_surfaceome(self):
        
        self.receptors_by_resource['Surfaceome'] = (
            self.annot.annots['Surfaceome'].get_subset(
                mainclass = 'Receptors',
            )
        )
    
    
    def add_ecm_matrisome(self):
        
        self.ecm_by_resource['Matrisome'] = (
            self.annot.annots['Matrisome'].get_subset(
                mainclass = 'Core matrisome',
            )
        )
    
    
    def add_ecm_go(self):
        
        self.ecm_by_resource['GO'] = (
            self.annot.annots['GO_Intercell'].get_subset(
                mainclass = 'ecm structure',
            )
        )
    
    
    def add_ecm_matrixdb(self):
        
        self.ecm_by_resource['MatrixDB'] = (
            self.annot.annots['MatrixDB_ECM'].to_set()
        )


class Intercell(session_mod.Logger):
    
    
    def __init__(self, network):
        
        session_mod.Logger.__init__(self, name = 'intercell')
        
        self.network = network
    
    
    def reload(self):
        """
        Reloads the object from the module level.
        """
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
