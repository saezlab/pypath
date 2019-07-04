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
        'protein_name',
        'transmembrane',
        'peripheral',
        'secreted',
        'secreted_desc',
        'secreted_highlight',
        'receptor',
        'receptor_desc',
        'integrin',
        'other',
        'other_desc',
        'tags',
        'tags_reason',
        'tags_description',
    ],
)


CellPhoneDBInteraction = collections.namedtuple(
    'CellPhoneDBInteraction',
    [
        'id_cp_interaction',
        'partner_a',
        'partner_b',
        'protein_name_a',
        'protein_name_b',
        'annotation_strategy',
        'source',
    ],
)


CellPhoneDBComplex = collections.namedtuple(
    'CellPhoneDBComplex',
    [
        'complex_name',
        'uniprot_1',
        'uniprot_2',
        'uniprot_3',
        'uniprot_4',
        'transmembrane',
        'peripheral',
        'secreted',
        'secreted_desc',
        'secreted_highlight',
        'receptor',
        'receptor_desc',
        'integrin',
        'other',
        'other_desc',
        'pdb_id',
        'pdb_structure',
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


categories_receiver = {
    'receptor',
    'transporter',
    'gap_junction',
    'tight_junction',
    'adhesion',
}


categories = categories_transmitter | categories_receiver


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
            os.path.join(self.output_dir, '%s_input.csv' % name),
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
        
        
        def get_id_name(entity):
            
            id_ = entity.__str__()
            
            name = (
                id_
                    if 'COMPLEX' in id_ else
                mapping.map_name0(
                    id_,
                    'uniprot',
                    'uniprot-entry',
                )
            )
            
            return id_, name
        
        
        iuphar = {'Guide2Pharma', 'Guide2Pharma_CP'}
        
        self._entities = set()
        self.interaction = set()
        
        int_id = 0
        for iaction in self.network.records.itertuples():
            
            if (
                # either the interaction is directed and the source is an
                # intercellular communication transmitter while the target
                # is a receiver
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
                categories_receiver
            ):
                
                id_a, name_a = get_id_name(iaction.id_a)
                id_b, name_b = get_id_name(iaction.id_b)
                
                self.interaction.add(
                    CellPhoneDBInteraction(
                        id_cp_interaction = 'CPI-%06u' % int_id,
                        partner_a = id_a,
                        partner_b = id_b,
                        protein_name_a = name_a,
                        protein_name_b = name_b,
                        annotation_strategy = (
                            'OmniPath,%s' % ','.join(sorted(iaction.sources))
                        ),
                        source = ','.join(
                            'PMID: %s' % pmid
                            for pmid in sorted(iaction.references)
                        ),
                    )
                )
                
                self._entities.add(iaction.id_a)
                self._entities.add(iaction.id_b)
                
                int_id += 1
    
    
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
                        protein_name = mapping.map_name0(
                            comp,
                            'uniprot',
                            'uniprot-entry',
                        ),
                        transmembrane = 'transmembrane' in classes,
                        peripheral = 'cell_surface' in classes,
                        secreted = 'secreted' in classes,
                        secreted_desc = '',
                        secreted_highlight = '',
                        receptor = 'receptor' in classes,
                        receptor_desc = '',
                        integrin = comp in integrins,
                        other = '',
                        other_desc = '',
                        tags = '',
                        tags_reason = '',
                        tags_description = '',
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
                    complex_name = entity.__str__(),
                    uniprot_1 = get_component(components, 0),
                    uniprot_2 = get_component(components, 1),
                    uniprot_3 = get_component(components, 2),
                    uniprot_4 = get_component(components, 3),
                    transmembrane = 'transmembrane' in classes,
                    peripheral = 'cell_surface' in classes,
                    secreted = 'secreted' in classes,
                    secreted_desc = '',
                    secreted_highlight = '',
                    receptor = 'receptor' in classes,
                    receptor_desc = '',
                    integrin = entity in integrins,
                    other = '',
                    other_desc = '',
                    pdb_id = (
                        ';'.join(sorted(entity.ids['PDB']))
                            if 'PDB' in entity.ids else
                        ''
                    ),
                    pdb_structure = '',
                    stoichiometry = entity.stoichiometry_str_genesymbols,
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
                ensembl_genes = mapping.map_name(comp, 'uniprot', 'ensembl')
                
                for ensembl in ensembl_genes:
                    
                    self.gene.add(
                        CellPhoneDBGene(
                            gene_name = name,
                            uniprot = comp,
                            hgnc_symbol = name,
                            ensembl = ensembl,
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
