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
import pypath.annot as annot
import pypath.intercell_annot as intercell_annot


IntercellRole = collections.namedtuple(
    'IntercellRole',
    ['source', 'role'],
)



class IntercellAnnotation(annot.CustomAnnotation):
    
    
    def __init__(self, class_definitions = None):
        
        class_definitions = (
            class_definitions or intercell_annot.default_intercell_classes
        )
        
        annot.CustomAnnotation.__init__(
            self,
            class_definitions = class_definitions,
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
        """
        Creates a consensus annotation of extracellular matrix structural
        protein components.
        """
        
        self.ecm_by_resource = {}
        
        self.add_ecm_go()
        self.add_ecm_matrisome()
        self.add_ecm_matrixdb()
        
        self.ecm_proteins = set.union(*self.ecm_by_resource.values())
    
    
    def collect_ligands(self):
        """
        Creates a consensus annotation of protein ligands in intercellular
        communication.
        """
        
        self.ligands_by_resource = {}
        
        self.add_ligands_cellphonedb()
        self.add_ligands_go()
        
        self.ligands = set.union(*self.ligands_by_resource.values())
    
    
    def collect_ec_enzymes(self):
        """
        Creates a consensus annotation of extracellular enzymes.
        """
        
        self.ec_enzymes = {}
    
    
    def collect_extracellular(self):
        
        self.extracellular_by_resource = {}
        
        
        self.add_extracellular_locatome()
        self.add_extracellular_go()
        self.add_extracellular_matrixdb()
        
        self.extracellular = set.union(
            *self.extracellular_by_resource.values()
        )
    
    (
        
    )
    
    def add_receptors_hpmr(self):
        
        self.receptors_by_resource['HPMR'] = (
            annot.db.annots['HPMR'].to_set()
        )
    
    
    def add_receptors_cellphonedb(self):
        
        self.receptors_by_resource['CellPhoneDB'] = (
            annot.db.annots['CellPhoneDB'].get_subset(
                receptor = bool,
                transmembrane = True,
            )
        )
    
    
    def add_receptors_go(self):
        
        self.receptors_by_resource['GO'] = (
            annot.db.annots['GO_Intercell'].get_subset(
                mainclass = 'receptors',
            )
        )
    
    
    def add_receptors_surfaceome(self):
        
        self.receptors_by_resource['Surfaceome'] = (
            annot.db.annots['Surfaceome'].get_subset(
                mainclass = 'Receptors',
            )
        )
    
    
    def add_ecm_matrisome(self):
        
        self.ecm_by_resource['Matrisome'] = (
            annot.db.annots['Matrisome'].get_subset(
                mainclass = 'Core matrisome',
            )
        )
    
    
    def add_ecm_go(self):
        
        self.ecm_by_resource['GO'] = (
            annot.db.annots['GO_Intercell'].get_subset(
                mainclass = 'ecm structure',
            )
        )
    
    
    def add_ecm_matrixdb(self):
        
        self.ecm_by_resource['MatrixDB'] = (
            annot.db.annots['MatrixDB_ECM'].to_set()
        )
    
    
    def add_ligands_cellphonedb(self):
        
        self.ligands_by_resource['CellPhoneDB'] = (
            annot.db.annots['CellPhoneDB'].get_subset(
                secreted = bool,
            )
        )
    
    
    def add_ligands_go(self):
        
        self.ligands_by_resource['GO'] = (
            annot.db.annots['GO_Intercell'].get_subset(
                mainclass = 'ligands',
            )
        )
    
    
    def add_extracellular_locatome(self):
        
        self.extracellular_by_resource['Locatome'] = (
            annot.db['Locatome'].get_subset(
                location = {'extracellular', 'extracellular region'}
            )
        ) & (
            annot.db['Locatome'].get_subset(cls = 'secretome')
        )
    
    
    def add_extracellular_go(self):
        
        self.extracellular_by_resource['GO_Intercell'] = (
            annot.db['GO_Intercell'].get_subset(mainclass = 'extracellular')
        )
    
    
    def add_extracellular_matrixdb(self):
        
        self.extracellular_by_resource['MatrixDB'] = set.union(
            annot.db.annots['MatrixDB_Secreted'].to_set(),
            annot.db.annots['MatrixDB_ECM'].to_set(),
        )
    
    
    def add_extracellular_surfaceome(self):
        
        self.extracellular_by_resource['Surfaceome'] = (
            annot.db.annots['Surfaceome'].to_set()
        )
    
    
    def add_extracellular_membranome(self):

        self.extracellular_by_resource['Membranome'] = (
            annot.db.annots['Membranome'].to_set()
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
