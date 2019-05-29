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

import os
import imp
import collections

import pandas as pd

import pypath.mapping as mapping
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


categories_transmitter = {
    'ligand',
    'gap_junction',
    'tight_junction',
    'extracellular_enzyme',
}


categories_reciever = {
    'receptor',
    'transporter',
    'gap_junction',
    'tight_junction',
    'adhesion',
}


categories = categories_transmitter | categories_reciever


class CellPhoneDB(session_mod.Logger):
    
    def __init__(
            self,
            output_dir = None,
            network = None,
            annotation = None,
            network_param = None,
            annot_param = None,
            omnipath = False,
        ):
        
        if not hasattr(self, '_log_name'):
            
            session_mod.Logger.__init__(self, name = 'cellphonedb')
        
        self.output_dir = output_dir
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
        
        self.setup()
        self.build()
        self.export()
    
    
    def setup(self):
        
        self.setup_paths()
        self.setup_network()
        self.setup_annotation()
    
    
    def setup_paths(self):
        
        self.output_dir = self.output_dir or '.'
        self._set_path('gene')
        self._set_path('protein')
        self._set_path('complex')
        self._set_path('interaction')
    
    
    def _set_path(self, name):
        
        setattr(
            self,
            '%s_path' % name,
            os.path.join(self.output_dir, '%s.csv' % name),
        )
    
    
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
    
    
    def build(self):
        
        self.build_interaction()
        self.build_protein()
        self.build_complex()
        self.build_gene()
        self.create_dataframes()
    
    
    def build_interaction(self):
        
        iuphar = {'Guide2Pharma', 'Guide2Pharma_CP'}
        
        self._entities = set()
        self.interaction = set()
        
        for iaction in self.network.records.itertuples():
            
            if (
                # either the interaction is directed and the source is an
                # intercellular communication transmitter while the target
                # is a reciever
                iaction.directed and
                (
                    self.annotation.classes_by_entity(iaction.id_a) &
                    categories_transmitter
                ) and
                (
                    self.annotation.classes_by_entity(iaction.id_b) &
                    categories_transmitter
                )
            ) or (
                # or it is indirected but both partners belong to the same
                # intercellular communication role category
                self.annotation.classes_by_entity(iaction.id_a) &
                self.annotation.classes_by_entity(iaction.id_b) &
                categories_transmitter
            ) or (
                self.annotation.classes_by_entity(iaction.id_a) &
                self.annotation.classes_by_entity(iaction.id_b) &
                categories_reciever
            ):
                
                self.interaction.add(
                    CellPhoneDBInteraction(
                        comments_interaction = '',
                        dlrp = False, # we don't have this info
                        family = '',
                        iuphar = bool(iaction.sources & iuphar),
                        multidata_name_1 = iaction.id_a.__str__(),
                        multidata_name_2 = iaction.id_b.__str__(),
                        score_1 = '',
                        score_2 = '',
                        source = ';'.join(sorted(iaction.sources)),
                    )
                )
                
                self._entities.add(iaction.id_a)
                self._entities.add(iaction.id_b)
    
    
    def build_protein(self):
        
        integrins = annot.db.annots['Integrins']
        
        self.protein = set()
        
        for entity in self._entities:
            
            # we add the components of the complexes to the protein data
            # frame; I don't know if it's necessary but does not harm I guess
            if hasattr(entity, 'components'):
                
                components = entity.components
                
            else:
                
                components = (entity,)
            
            for comp in components:
                
                classes = self.annotation.classes_by_entity(comp)
                
                self.protein.add(
                    CellPhoneDBProtein(
                        uniprot = comp.__str__(),
                        adhesion = 'adhesion' in classes,
                        cytoplasm = 'intracellular' in classes,
                        entry_name = mapping.map_name0(
                            comp,
                            'uniprot',
                            'uniprot-entry',
                        ),
                        extracellular = 'extracellular' in classes,
                        integrin_interaction = comp in integrins,
                        other = '',
                        other_desc = '',
                        pdb_id = '',
                        pdb_structure = '',
                        peripheral = 'cell_surface' in classes,
                        receptor = 'receptor' in classes,
                        receptor_desc = '',
                        secreted_desc = '',
                        secreted_highlight = '',
                        secretion = 'secreted' in classes,
                        stoichiometry = '',
                        tags_description = '',
                        transmembrane = 'transmembrane' in classes,
                        transporter = 'transporter' in classes,
                        tags = '',
                        tags_reason = '',
                    )
                )
    
    
    def build_complex(self):
        
        def get_component(components, idx):
            
            return components[idx] if len(components) > idx else ''
        
        integrins = annot.db.annots['Integrins']
        
        self.complex = set()
        
        for entity in self._entities:
            
            if not hasattr(entity, 'components'):
                
                continue
            
            # hey, what to do with complexes with more than 4 components???
            components = sorted(entity.components.keys())
            
            classes = self.annotation.classes_by_entity(entity)
            
            self.complex.add(
                CellPhoneDBComplex(
                    name = entity.__str__(),
                    uniprot_1 = get_component(components, 0),
                    uniprot_2 = get_component(components, 1),
                    uniprot_3 = get_component(components, 2),
                    uniprot_4 = get_component(components, 3),
                    adhesion = 'adhesion' in classes,
                    cytoplasm = 'intracellular' in classes,
                    extracellular = 'extracellular' in classes,
                    integrin_interaction = entity in integrins,
                    other = '',
                    other_desc = '',
                    peripheral = 'cell_surface' in classes,
                    receptor = 'receptor' in classes,
                    receptor_desc = '',
                    secreted_desc = '',
                    secreted_highlight = '',
                    secretion = 'secreted' in classes,
                    transmembrane = 'transmembrane' in classes,
                    transporter = 'transporter' in classes,
                    pdb_structure = '',
                    pdb_id = (
                        ';'.join(sorted(entity.ids['PDB']))
                            if 'PDB' in entity.ids else
                        ''
                    ),
                    stoichiometry = entity.stoichiometry,
                    comments_complex = '',
                )
            )
    
    
    def build_gene(self):
        
        self.gene = set()
        
        for entity in self._entities:
            
            # we add the components of the complexes to the protein data
            # frame; I don't know if it's necessary but does not harm I guess
            if hasattr(entity, 'components'):
                
                components = entity.components
                
            else:
                
                components = (entity,)
            
            for comp in components:
                
                name = mapping.map_name0(comp, 'uniprot', 'genesymbol')
                
                self.gene.add(
                    CellPhoneDBGene(
                        gene_name = name,
                        uniprot = comp,
                        hgnc_symbol = name,
                        ensembl = '',
                    )
                )
    
    
    def create_dataframes(self):
        
        for name in ('interaction', 'protein', 'complex', 'gene'):
            
            data = list(getattr(self, name))
            columns = globals()['CellPhoneDB%s' % name.capitalize()]._fields
            
            setattr(
                self,
                '%s_dataframe' % name,
                pd.DataFrame(
                    data,
                    columns = columns,
                )
            )
    
    
    def export(self):
        
        self.export_gene()
        self.export_protein()
        self.export_complex()
        self.export_interaction()
    
    
    def _export(self, name):
        
        path = getattr(self, '%s_path' % name)
        getattr(self, '%s_dataframe' % name).to_csv(
            path,
            sep = ',',
            index = False,
        )
        self._log('Data frame `%s` exported to `%s`.' % (name, path))
    
    
    def export_gene(self):
        
        self._export('gene')
    
    
    def export_protein(self):
        
        self._export('protein')
    
    
    def export_complex(self):
        
        self._export('complex')
    
    
    def export_interaction(self):
        
        self._export('interaction')
