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

from future.utils import iteritems

import os
import importlib as imp
import collections
import copy

import pandas as pd

import pypath.utils.mapping as mapping
import pypath.share.session as session_mod
import pypath.legacy.main as main_mod
import pypath.core.network as network_mod
import pypath.resources.data_formats as data_formats
import pypath.core.annot as annot
import pypath.core.intercell as intercell
import pypath.omnipath.legacy.legacy as omnipath


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
            complexes = None,
            intercell = None,
            network_param = None,
            intercell_param = None,
            use_omnipath = False,
            network_pickle = None,
            annotation_pickle = None,
            intercell_pickle = None,
            complex_pickle = None,
        ):
        
        session_mod.Logger.__init__(self, name = 'cellphonedb')
        
        self.output_dir = output_dir
        self.network = network
        self.annot = annotation
        self.complex = complexes
        self.intercell = intercell
        self.network_param = network_param or {}
        self.intercell_param = intercell_param or {}
        self.use_omnipath = use_omnipath
        
        omnipath.OmniPath.__init__(
            self,
            network_pickle = network_pickle,
            annotation_pickle = annotation_pickle,
            intercell_pickle = intercell_pickle,
            complex_pickle = complex_pickle,
            load_network = False, # we load the network here below
            load_complexes = not self.complex,
            load_annotations = not self.intercell and not self.annot,
            load_intercell = not self.intercell,
            load_enz_sub = False,
        )
    
    
    def reload(self):
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        import importlib as imp
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    
    def main(self):
        
        omnipath.OmniPath.main(self)
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
            
            if not self.use_omnipath and 'lst' not in self.network_param:
                
                resources = copy.deepcopy(data_formats.pathway_all)
                resources.update(copy.deepcopy(data_formats.ligand_receptor))
                self.network_param['lst'] = resources
            
            self.network = main_mod.PyPath()
            self.network.init_network(**self.network_param)
        
        if isinstance(self.network, main_mod.PyPath):
            
            self.network = network_mod.Network(records = self.network)
            
        if isinstance(self.network, network_mod.Network):
            
            pass
    
    
    def setup_annotation(self):
        
        if self.intercell is None:
            
            self.intercell = intercell.IntercellAnnotation(
                **self.intercell_param
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
        self.cpdb_interaction = {}
        
        int_id = 0
        for iaction in self.network.records.itertuples():
            
            if (
                # either the interaction is directed and the source is an
                # intercellular communication transmitter while the target
                # is a receiver
                iaction.directed and
                (
                    self.intercell.classes_by_entity(iaction.id_a) &
                    categories_transmitter
                ) and
                (
                    self.intercell.classes_by_entity(iaction.id_b) &
                    categories_transmitter
                )
            ) or (
                # or it is undirected but both partners belong to the same
                # intercellular communication role category
                self.intercell.classes_by_entity(iaction.id_a) &
                self.intercell.classes_by_entity(iaction.id_b) &
                categories_transmitter
            ) or (
                self.intercell.classes_by_entity(iaction.id_a) &
                self.intercell.classes_by_entity(iaction.id_b) &
                categories_receiver
            ):
                
                id_a, name_a = get_id_name(iaction.id_a)
                id_b, name_b = get_id_name(iaction.id_b)
                
                key = (id_a, id_b)
                
                int_id_carry_over = (
                    int(self.cpdb_interaction[key].id_cp_interaction[5:])
                        if key in self.cpdb_interaction else
                    None
                )
                
                resources_carry_over, refs_carry_over = (
                    (
                        set(
                            self.cpdb_interaction[
                                key
                            ].annotation_strategy.split(',')
                        ),
                        set(
                            pmid[7:]
                            for pmid in
                            self.cpdb_interaction[key].source.split(',')
                        )
                    )
                    if key in self.cpdb_interaction else
                    (set(), set())
                )
                
                self.cpdb_interaction[key] = (
                    CellPhoneDBInteraction(
                        id_cp_interaction = 'CPI-%06u' % (
                            int_id_carry_over or int_id
                        ),
                        partner_a = id_a,
                        partner_b = id_b,
                        protein_name_a = name_a,
                        protein_name_b = name_b,
                        annotation_strategy = (
                            ','.join(sorted(
                                resources_carry_over |
                                iaction.sources |
                                {'OmniPath'}
                            ))
                        ),
                        source = ','.join(
                            'PMID: %s' % pmid
                            for pmid in sorted(
                                refs_carry_over |
                                (iaction.references or set())
                            )
                        ),
                    )
                )
                
                self._entities.add(iaction.id_a)
                self._entities.add(iaction.id_b)
                
                if int_id_carry_over is None:
                    
                    int_id += 1
        
        self.cpdb_interaction = sorted(
            self.cpdb_interaction.values(),
            key = lambda i: i[0]
        )
    
    
    def build_protein(self):
        
        integrins = annot.db.annots['Integrins']
        
        self.cpdb_protein = set()
        
        for entity in self._entities:
            
            # we add the components of the complexes to the protein data
            # frame; I don't know if it's necessary but does not harm I guess
            if hasattr(entity, 'components'):
                
                components = entity.components
                
            else:
                
                components = (entity,)
            
            for comp in components:
                
                classes = self.intercell.classes_by_entity(comp)
                
                self.cpdb_protein.add(
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
        
        self.cpdb_complex = set()
        
        for entity in self._entities:
            
            if not hasattr(entity, 'components'):
                
                continue
            
            # hey, what to do with complexes with more than 4 components???
            components = sorted(entity.components.keys())
            
            classes = self.intercell.classes_by_entity(entity)
            
            self.cpdb_complex.add(
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
        
        self.cpdb_gene = set()
        
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
                    
                    self.cpdb_gene.add(
                        CellPhoneDBGene(
                            gene_name = name,
                            uniprot = comp,
                            hgnc_symbol = name,
                            ensembl = ensembl,
                        )
                    )
    
    
    def create_dataframes(self):
        
        for name in ('interaction', 'protein', 'complex', 'gene'):
            
            self._log('Creating CellPhoneDB data frame `%s`.' % name)
            
            data = list(getattr(self, 'cpdb_%s' % name))
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


# example:
# IL12 receptor: 'COMPLEX:P42701-Q99665'
# IL12 ligand: 'COMPLEX:P29459-P29460'
# 
