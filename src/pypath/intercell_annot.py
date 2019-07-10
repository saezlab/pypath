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

import pypath.go as go
import pypath.dataio as dataio
import pypath.annot_formats as af


"""
Gene Ontology annotations to select categories relevant in intercellular
signaling.
"""
go_combined_classes = {
    'junction':
        """
        cell junction OR
        cell junction assembly OR
        cell junction organization OR
        intercellular bridge organization OR
        gap junction-mediated intercellular transport OR
        gap junction
        """,
    'gap junction':
        """
        gap junction
        """,
    'tight junction':
        """
        tight junction
        """,
    'cell-substrate junction':
        """
        cell-substrate junction
        """,
    'cell-cell junction':
        """
        cell-cell junction
        """,
    'extracellular':
        """
        extracellular region OR
        extracellular region part
        """,
    'intracellular':
        """
        intracellular organelle OR
        intracellular organelle lumen OR
        intracellular
        """,
    'cell_surface':
        """
        cell surface OR
        external side of plasma membrane
        """,
    'transmembrane':
        """
        integral to membrane OR
        transmembrane signaling receptor activity
        """,
    'ecm':
        """
        extracellular matrix OR
        complex of collagen trimers OR
        collagen network OR
        banded collagen fibril OR
        collagen beaded filament OR
        elastic fiber OR
        fibronectin fibril
        """,
    'extracell enzyme':
        """
        catalytic activity AND
        extracellular region
        """,
    'enzyme':
        """
        catalytic activity
        """,
    'extracell peptidase':
        """
        peptidase activity AND
        extracellular region
        """,
    'peptidase':
        """
        peptidase activity
        """,
    'growth factor binding':
        """
        growth factor binding AND
        extracellular space
        """,
    'receptor regulation':
        """
        (receptor regulator activity OR
        regulation of receptor recycling OR
        receptor clustering OR
        receptor diffusion trapping OR
        membrane raft localization) AND
        (extracellular region OR
        cell surface OR
        external side of plasma membrane OR
        intrinsic component of plasma membrane)
        """,
    'receptor inhibition':
        """
        receptor inhibitor activity AND
        (extracellular region OR
        cell surface OR
        external side of plasma membrane OR
        intrinsic component of plasma membrane)
        """,
    'receptor activation':
        """
        signaling receptor activator activity AND
        (extracellular region OR
        cell surface OR
        external side of plasma membrane OR
        intrinsic component of plasma membrane)
        """,
    'ligands':
        """
        receptor ligand activity AND
        (extracellular region OR
        cell surface OR
        external side of plasma membrane)
        """,
    'secreted ligands':
        """
        receptor ligand activity AND
        extracellular region
        """,
    'surface ligands':
        """
        receptor ligand activity AND
        (cell surface OR
        external side of plasma membrane)
        """,
    'receptors':
        """
        (signaling receptor activity OR
        cell surface receptor signaling pathway OR
        transmembrane signaling receptor activity OR
        receptor complex OR
        (cellular response to stimulus AND signal transduction)) AND
        (cell surface OR
        external side of plasma membrane)
        """,
    'membrane ligands':
        """
        receptor ligand activity AND
        (cell surface OR
        external side of plasma membrane)
        """,
    'hormone receptors':
        """
        hormone binding AND
        (cell surface OR
        external side of plasma membrane)
        """,
    'ecm structure':
        """
        (extracellular region OR
        extracellular matrix) AND
        (structural molecule activity OR
        extracellular matrix structural constituent OR
        structural constituent of bone)
        """,
    'ecm production':
        """
        extracellular matrix assembly OR
        extracellular matrix organization OR
        extracellular matrix constituent secretion OR
        collagen metabolic process
        """,
    'endocytosis':
        """
        cargo adaptor activity
        """,
    'adhesion to matrix':
        """
        cell adhesion mediator activity OR
        extracellular matrix binding OR
        hydroxyapatite binding OR
        cell-substrate adhesion
        """,
    'response to adhesion':
        """
        cellular response to cell-matrix adhesion OR
        contact inhibition OR
        establishment or maintenance of cell polarity OR
        extracellular matrix-cell signaling
        """,
    'adhesion to other cells':
        """
        cell-cell adhesion
        """,
    'adhesion':
        """
        cell adhesion mediator activity OR
        extracellular matrix binding OR
        hydroxyapatite binding OR
        cell adhesion OR
        cell-cell adhesion OR
        cell-substrate adhesion OR
        cell adhesion molecule production OR
        cell-cell adhesion in response to extracellular stimulus OR
        cellular response to cell-matrix adhesion OR
        contact inhibition OR
        establishment or maintenance of cell polarity
        """,
    'cell-cell signaling':
        """
        cell-cell signaling
        """,
    'transport':
        """
        intrinsic component of plasma membrane AND
        transmembrane transporter activity
        """,
    'regulation of transport':
        """
        (cell surface OR
        external side of plasma membrane OR
        extracellular region) AND
        (regulation of transmembrane transporter activity OR
        channel regulator activity)
        """,
    'ion channels':
        """
        (ion channel activity OR
        ion transmembrane transporter activity) AND
        intrinsic component of plasma membrane
        """,
    'regulation of ion channels':
        """
        (cell surface OR
        external side of plasma membrane OR
        extracellular region) AND
        (ion channel regulator activity OR
        regulation of ion transmembrane transporter activity OR
        channel regulator activity)
        """,
    'autocrine signaling':
        """
        autocrine signaling
        """,
    'paracrine signaling':
        """
        paracrine signaling
        """,
    'signal release':
        """
        signal release OR
        hormone secretion OR
        hormone metabolic process OR
        cytokine production OR
        cytokine secretion
        """,
    'secretion':
        """
        exocytic process OR
        secretion by cell
        """,
    'communication by exosomes':
        """
        extracellular vesicle biogenesis OR
        extracellular exosome assembly OR
        vesicle-mediated intercellular transport OR
        cell-cell signaling via exosome
        """,
    'membrane docking':
        """
        membrane docking AND
        (extracellular region OR
        cell surface OR
        external side of plasma membrane)
        """
}


annotation_categories = [
    ('MatrixDB_Secreted',),
    ('MatrixDB_ECM',),
    ('Surfaceome',),
    ('GO_Intercell', 'adhesion'),
    ('GO_Intercell', 'cell_surface'),
    ('GO_Intercell', 'ecm'),
    ('GO_Intercell', 'ecm structure'),
    ('GO_Intercell', 'extracell enzyme'),
    ('GO_Intercell', 'extracell peptidase'),
    ('GO_Intercell', 'extracellular'),
    ('GO_Intercell', 'hormone receptors'),
    ('GO_Intercell', 'ion channels'),
    ('GO_Intercell', 'junction'),
    ('GO_Intercell', 'ligands'),
    ('GO_Intercell', 'receptors'),
    ('MatrixDB_Membrane',),
    ('Membranome', 'Plasma membrane', 'extracellular side'),
    ('CSPA',),
    ('Locate', 'extracellular'),
    ('Locate', 'extracellular region'),
    ('Locate', 'plasma membrane'),
    ('Matrisome', 'Core matrisome'),
    ('Matrisome', 'Matrisome-associated'),
    ('Matrisome', 'Core matrisome', 'Collagens'),
    ('Matrisome', 'Core matrisome', 'ECM Glycoproteins'),
    ('Matrisome', 'Core matrisome', 'Proteoglycans'),
    ('Matrisome', 'Matrisome-associated', 'ECM Regulators'),
    ('Matrisome', 'Matrisome-associated', 'ECM-affiliated Proteins'),
    ('Matrisome', 'Matrisome-associated', 'Secreted Factors'),
]

[
    ('GO_Intercell', 'adhesion'),
    ('GO_Intercell', 'cell_surface'),
    ('GO_Intercell', 'ecm structure'),
    ('GO_Intercell', 'extracell enzyme'),
    ('GO_Intercell', 'extracell peptidase'),
    ('GO_Intercell', 'ion channels'),
    ('GO_Intercell', 'junction'),
    ('GO_Intercell', 'ligands'),
    ('GO_Intercell', 'receptors'),
]


go_single_terms = {
    
    # cellular component
    'C': {
        # junction
        'cell junction',
        
        # extracellular
        'extracellular region',
        'extracellular region part',
        
        # extracellular_matrix
        'extracellular matrix',
        'complex of collagen trimers',
        'collagen network',
        'banded collagen fibril',
        'collagen beaded filament',
        'elastic fiber',
        'fibronectin fibril',
        
        # exosome
        # could not find sub-term for their membrane or lumen
        'extracellular vesicle',
        
        # cell_surface
        # only plasma membrane components facing outside
        'cell surface',
        'external side of plasma membrane',
        
        # these contains also intracellular components
        'plasma membrane',
        'extrinsic component of plasma membrane',
        'intrinsic component of plasma membrane',
        'cytoplasmic side of plasma membrane',
        
        'immunological synapse',
        
        'clathrin-coated pit',
        'plasma membrane raft',
        
        'presynaptic membrane',
        'presynaptic endocytic zone',
        'presynaptic endocytic zone membrane',
        'intrinsic component of presynaptic membrane',
        'extrinsic component of presynaptic membrane',
        'presynaptic active zone membrane',
        
        'postsynaptic density membrane',
        'intrinsic component of postsynaptic density membrane',
        'extrinsic component of postsynaptic density membrane',
        
        'synaptic vesicle',
        'intrinsic component of synaptic vesicle membrane',
        'synaptic vesicle lumen',
        'synaptic vesicle membrane',
        'extracellular matrix of synaptic cleft',
        
        'axolemma',
        'neuron projection membrane',
        'neuronal cell body membrane',
        'photoreceptor inner segment membrane',
        'photoreceptor outer segment membrane',
        'stereocilium membrane',
        'stereocilia coupling link',
        
    },
    
    # molecular function
    'F': {
        
        # upper level terms which help to categorize
        # molecular functions in the intercellular signaling
        'molecular carrier activity',
        'cargo receptor activity',
        'binding',
        'regulation of binding',
        'positive regulation of binding',
        'negative regulation of binding',
        'protein folding chaperone',
        'antioxidant activity',
        
        # upper level generic terms for regulation
        'regulation of molecular function',
        'negative regulation of molecular function',
        'positive regulation of molecular function',
        
        # enzymes
        'catalytic activity',
        'regulation of catalytic activity',
        'positive regulation of catalytic activity',
        'negative regulation of catalytic activity',
        'enzyme regulator activity',
        'enzyme activator activity',
        'enzyme inhibitor activity',
        'catalytic activity, acting on a protein',
        # peptidases
        'peptidase activity',
        'peptidase regulator activity',
        'peptidase activator activity',
        'peptidase inhibitor activity',
        'regulation of peptidase activity',
        'peptidase activator activity',
        # receptors
        'receptor regulator activity',
        'signaling receptor activator activity',
        'receptor inhibitor activity',
        'receptor ligand activity',
        'neurotransmitter receptor regulator activity',
        'signaling receptor activity',
        'negative regulation of signaling receptor activity',
        'positive regulation of signaling receptor activity',
        'signaling receptor activator activity',
        'receptor inhibitor activity',
        'regulation of signaling receptor activity',
        'receptor complex',
        (
            'transforming growth factor beta receptor,'
            'cytoplasmic mediator activity'
        ),
        # other relevant binding activities
        'antigen binding',
        'hormone binding',
        'neurotransmitter binding',
        
        # endocytosis
        'cargo adaptor activity',
        
        # ECM, structural proteins
        'structural molecule activity',
        'extracellular matrix structural constituent',
        'structural constituent of bone',
        
        # adhesion to base membrane and ECM
        'cell adhesion mediator activity',
        'extracellular matrix binding',
        'hydroxyapatite binding',
        
        # transporters
        'transporter activity',
        'regulation of transporter activity',
        'positive regulation of transporter activity',
        'negative regulation of transporter activity',
        'transmembrane transporter activity',
        'drug transmembrane transporter activity',
        'regulation of transmembrane transporter activity',
        # channels
        'channel regulator activity',
        'channel activator activity',
        'channel inhibitor activity',
        # ion channels
        'ion channel inhibitor activity',
        'ion channel regulator activity',
        'ion transmembrane transporter activity',
        'regulation of ion transmembrane transporter activity',
        'positive regulation of ion transmembrane transporter activity',
        'negative regulation of ion transmembrane transporter activity',
    },
    
    # biological process
    'P': {
        
        # cell adhesion
        'cell adhesion',
        'cell-cell adhesion',
        'cell-substrate adhesion',
        'cell adhesion molecule production',
        'cell-cell adhesion in response to extracellular stimulus',
        'cellular response to cell-matrix adhesion',
        'contact inhibition',
        'establishment or maintenance of cell polarity'
        
        # cellular responses
        'myofibroblast cell apoptotic process',
        'fibroblast apoptotic process',
        'epithelial cell apoptotic process',
        
        # cell activation
        'cell activation',
        'endothelial cell activation',
        'fibroblast activation',
        'leukocyte activation',
        
        # cell communication
        'cell communication',
        'cell communication by chemical coupling',
        'cell communication by electrical coupling',
        'cell-cell signaling',
        'autocrine signaling',
        'paracrine signaling',
        'cell-cell signaling via exosome',
        'epithelial-mesenchymal cell signaling',
        'cellular response to extracellular stimulus',
        'synaptic signaling',
        'cell-cell recognition',
        
        # signal release
        'signal release',
        'hormone secretion',
        'exocytic process',
        'secretion by cell',
        'hormone metabolic process',
        'cytokine production',
        'cytokine secretion',
        
        # receptors
        'regulation of receptor recycling',
        'receptor clustering',
        'receptor diffusion trapping',
        'membrane raft localization',
        
        # exosomes
        'extracellular vesicle biogenesis',
        'extracellular exosome assembly',
        'vesicle-mediated intercellular transport',
        
        # junction
        'cell junction assembly',
        'cell junction organization',
        'intercellular bridge organization',
        'gap junction-mediated intercellular transport',
        
        # ECM
        'extracellular matrix assembly',
        'extracellular matrix organization',
        'extracellular matrix constituent secretion',
        'extracellular matrix-cell signaling',
        'cell-matrix recognition',
        'collagen metabolic process',
        
        # motility
        'fibroblast migration',
        'epithelial cell migration',
        'endothelial cell migration',
        'leukocyte migration',
        'substrate-dependent cell migration',
        'cell chemotaxis to fibroblast growth factor',
        'endothelial cell chemotaxis',
        'fibroblast chemotaxis',
        'leukocyte chemotaxis',
        'epithelial structure maintenance',
        'connective tissue replacement', # this is fibrosis actually
        
        # channels
        'ion channel activity',
        
        # docking
        'membrane docking',
        'protein to membrane docking',
        'membrane to membrane docking',
        
        # kidney specific stuff
        'nephron',
        'outer medulla of kidney',
        'inner medulla of kidney',
        'kidney pyramid',
        
    },
    
}


"""
Higher level classes of intercellular communication roles.
"""
annot_combined_classes = (
    # receptor
    af.AnnotDef(
        name = 'receptor',
        source = af.AnnotOp(
            annots = (
                'receptor_cellphonedb',
                'receptor_surfaceome',
                'receptor_go',
                'receptor_hpmr',
                'receptor_ramilowski',
                'receptor_kirouac',
                'receptor_guide2pharma',
                'receptor_hgnc',
            ),
            op = set.union,
        ),
    ),
    af.AnnotDef(
        name = 'receptor_hgnc',
        source = af.AnnotOp(
            annots = (
                'interleukin_receptors_hgnc',
            ),
            op = set.union,
        ),
    ),
    af.AnnotDef(
        name = 'interleukin_receptors_hgnc',
        source = 'HGNC',
        args = {
            'mainclass': 'Interleukin receptors',
        },
    ),
    af.AnnotDef(
        name = 'receptor_cellphonedb',
        source = 'CellPhoneDB',
        args = {
            'receptor': bool,
            'transmembrane': True,
        },
    ),
    af.AnnotDef(
        name = 'receptor_go',
        source = 'GO_Intercell',
        args = {
            'mainclass': 'receptors',
        },
    ),
    af.AnnotDef(
        name = 'receptor_hpmr',
        source = 'HPMR',
        args = {
            'role': 'Receptor',
        },
    ),
    af.AnnotDef(
        name = 'receptor_surfaceome',
        source = 'Surfaceome',
        args = {
            'mainclass': 'Receptors',
        },
    ),
    af.AnnotDef(
        name = 'receptor_ramilowski',
        source = 'Ramilowski2015',
        args = {
            'mainclass': 'receptor',
        },
    ),
    af.AnnotDef(
        name = 'receptor_kirouac',
        source = 'Kirouac2010',
        args = {
            'mainclass': 'receptor',
        },
    ),
    af.AnnotDef(
        name = 'receptor_guide2pharma',
        source = 'Guide2Pharma',
        args = {
            'mainclass': 'receptor',
        },
    ),
    # ECM
    af.AnnotDef(
        name = 'ecm',
        source = af.AnnotOp(
            annots = (
                'ecm_matrixdb',
                'ecm_matrisome',
                'ecm_go',
            ),
            op = set.union,
        ),
    ),
    af.AnnotDef(
        name = 'ecm_matrisome',
        source = af.AnnotOp(
            annots = (
                af.AnnotDef(
                    name = 'ecm_matrisome_core',
                    source = 'Matrisome',
                    args = {
                        'mainclass': 'Core matrisome',
                    },
                ),
                af.AnnotOp(
                    annots = (
                        af.AnnotDef(
                            name = 'ecm_matrisome_affiliated',
                            source = 'Matrisome',
                            args = {
                                'mainclass': 'Matrisome-associated',
                                'subclass': 'ECM-affiliated Proteins',
                            },
                        ),
                        'cell_surface',
                    ),
                    op = set.difference,
                ),
            ),
            op = set.union,
        ),
    ),
    af.AnnotDef(
        name = 'ecm_matrixdb',
        source = 'MatrixDB',
        args = {
            'mainclass': 'ecm',
        },
    ),
    af.AnnotDef(
        name = 'ecm_go',
        source = 'GO_Intercell',
        args = {
            'mainclass': 'ecm structure',
        },
    ),
    # ligand
    af.AnnotDef(
        name = 'ligand',
        source = af.AnnotOp(
            annots = (
                'ligand_cellphonedb',
                'ligand_go',
                'ligand_hpmr',
                'ligand_ramilowski',
                'ligand_kirouac',
                'ligand_guide2pharma',
                'ligand_hgnc',
            ),
            op = set.union,
        ),
    ),
    af.AnnotDef(
        name = 'ligand_hgnc',
        source = af.AnnotOp(
            annots = (
                'interleukins_hgnc',
                'endogenous_ligands_hgnc',
                'chemokine_ligands_hgnc',
            ),
            op = set.union,
        ),
    ),
    af.AnnotDef(
        name = 'interleukins_hgnc',
        source = 'HGNC',
        args = {
            'mainclass': 'Interleukins',
        },
    ),
    af.AnnotDef(
        name = 'endogenous_ligands_hgnc',
        source = 'HGNC',
        args = {
            'mainclass': {
                'Endogenous ligands',
                'Neurotrophins',
                'GDNF family ligands',
            },
        },
    ),
    af.AnnotDef(
        name = 'chemokine_ligands_hgnc',
        source = 'HGNC',
        args = {
            'mainclass': 'Chemokine ligands',
        },
    ),
    af.AnnotDef(
        name = 'ligand_cellphonedb',
        source = 'CellPhoneDB',
        args = {
            'secreted': bool,
        },
    ),
    af.AnnotDef(
        name = 'ligand_go',
        source = 'GO_Intercell',
        args = {
            'mainclass': 'ligands',
        },
    ),
    af.AnnotDef(
        name = 'ligand_hpmr',
        source = 'HPMR',
        args = {
            'role': 'Ligand',
        },
    ),
    af.AnnotDef(
        name = 'ligand_ramilowski',
        source = 'Ramilowski2015',
        args = {
            'mainclass': 'ligand',
        },
    ),
    af.AnnotDef(
        name = 'ligand_kirouac',
        source = 'Kirouac2010',
        args = {
            'mainclass': 'ligand',
        },
    ),
    af.AnnotDef(
        name = 'ligand_guide2pharma',
        source = 'Guide2Pharma',
        args = {
            'mainclass': 'ligand',
        },
    ),
    # intracellular
    af.AnnotDef(
        name = 'intracellular',
        source = af.AnnotOp(
            annots = (
                'intracellular_locate',
                'intracellular_cellphonedb',
                'intracellular_comppi',
                'intracellular_go',
            ),
            op = set.union,
        ),
    ),
    af.AnnotDef(
        name = 'intracellular_locate',
        source = af.AnnotOp(
            annots = (
                af.AnnotDef(
                    name = 'locate_intracellular',
                    source = 'Locate',
                    args = {
                        'location': {
                            'centrosome',
                            'cytoplasm',
                            'endosomes',
                            'lysosomes',
                            'nucleus',
                            'plasma membrane',
                            'cytoplasmic membrane-bound vesicle',
                            'cytoplasmic vesicles',
                            'cytoskeleton',
                            'early endosomes',
                            'endoplasmic reticulum',
                            'golgi apparatus',
                            'er-golgi intermediate compartment',
                            'ergic',
                            'golgi cis cisterna',
                            'golgi medial cisterna',
                            'golgi trans cisterna',
                            'golgi trans face',
                            'inner mitochondrial membrane',
                            'late endosomes',
                            'lipid particles',
                            'medial-golgi',
                            'melanosome',
                            'microtubule',
                            'microtubule organizing center ',
                            'mitochondria',
                            'mitochondrial inner membrane',
                            'mitochondrial outer membrane',
                            'mitochondrion',
                            'nuclear envelope',
                            'nucleolus',
                            'nuclear speck',
                            'outer mitochondrial membrane',
                            'peroxisome',
                            'peroxisomes',
                            'sarcolemma',
                            'transport vesicle',
                        },
                    },
                ),
                af.AnnotDef(
                    name = 'locate_cytoplasmic',
                    source = 'Locate',
                    args = {
                        'cls': 'cytoplasmic',
                    },
                ),
            ),
            op = set.union,
        ),
    ),
    af.AnnotDef(
        name = 'intracellular_cellphonedb',
        source = 'CellPhoneDB',
        args = {
            'cytoplasm': bool,
        },
    ),
    af.AnnotDef(
        name = 'intracellular_comppi',
        source = 'ComPPI',
        args = {
            'location': {
                'cytosol',
                'nucleus',
                'mitochondrion',
            },
        },
    ),
    af.AnnotDef(
        name = 'intracellular_go',
        source = 'GO_Intercell',
        args = {
            'mainclass': 'intracellular',
        },
    ),
    # extracellular
    af.AnnotDef(
        name = 'extracellular',
        source = af.AnnotOp(
            annots = (
                'extracellular_locate',
                'extracellular_surfaceome',
                'extracellular_matrixdb',
                'extracellular_membranome',
                'extracellular_cspa',
                'extracellular_hpmr',
                'extracellular_cellphonedb',
                'extracellular_comppi',
            ),
            op = set.union,
        ),
    ),
    af.AnnotDef(
        name = 'extracellular_locate',
        source = af.AnnotOp(
            annots = (
                af.AnnotDef(
                    name = 'locate_extracellular',
                    source = 'Locate',
                    args = {
                        'location': {
                            'extracellular',
                            'extracellular region',
                        },
                    },
                ),
                af.AnnotDef(
                    name = 'locate_secretome',
                    source = 'Locate',
                    args = {
                        'cls': 'secretome',
                    },
                ),
            ),
            op = set.union,
        ),
    ),
    af.AnnotDef(
        name = 'extracellular_comppi',
        source = 'ComPPI',
        args = {
            'location': 'extracellular',
        },
    ),
    af.AnnotDef(
        name = 'extracellular_surfaceome',
        source = 'Surfaceome',
    ),
    af.AnnotDef(
        name = 'extracellular_matrixdb',
        source = 'Matrixdb',
    ),
    af.AnnotDef(
        name = 'extracellular_cspa',
        source = 'CSPA',
    ),
    af.AnnotDef(
        name = 'extracellular_membranome',
        source = 'Membranome',
        args = {
            'membrane': 'Plasma membrane',
            'side': 'extracellular side',
        },
    ),
    af.AnnotDef(
        name = 'extracellular_cellphonedb',
        source = 'CellPhoneDB',
        args = {
            'extracellular': bool,
        },
    ),
    af.AnnotDef(
        name = 'extracellular_hpmr',
        source = 'HPMR',
    ),
    af.AnnotDef(
        name = 'extracellular_membranome',
        source = 'Membranome',
        args = {
            'membrane': 'Plasma membrane',
            'side': 'extracellular side',
        },
    ),
    # cell surface
    af.AnnotDef(
        name = 'cell_surface',
        source = af.AnnotOp(
            annots = (
                'cell_surface_surfaceome',
                'cell_surface_go',
                'cell_surface_hpmr',
                'cell_surface_membranome',
                'cell_surface_cspa',
                'cell_surface_cellphonedb',
            ),
            op = set.union,
        ),
    ),
    af.AnnotDef(
        name = 'cell_surface_hpmr',
        source = 'HPMR',
        args = {
            'role': 'Receptor',
        },
    ),
    af.AnnotDef(
        name = 'cell_surface_surfaceome',
        source = 'Surfaceome',
    ),
    af.AnnotDef(
        name = 'cell_surface_cspa',
        source = 'CSPA',
    ),
    af.AnnotDef(
        name = 'cell_surface_cellphonedb',
        source = 'CellPhoneDB',
        args = {
            'method': lambda a: a.peripheral or a.transmembrane,
        },
    ),
    af.AnnotDef(
        name = 'cell_surface_go',
        source = 'GO_Intercell',
        args = {
            'mainclass': 'cell_surface',
        },
    ),
    af.AnnotDef(
        name = 'cell_surface_membranome',
        source = 'Membranome',
        args = {
            'membrane': 'Plasma membrane',
            'side': 'extracellular side',
        },
    ),
    # transmembrane
    af.AnnotDef(
        name = 'transmembrane',
        source = af.AnnotOp(
            annots = (
                'transmembrane_cellphonedb',
                'transmembrane_go',
                'transmembrane_opm',
                'transmembrane_locate',
                'transmembrane_topdb',
            ),
            op = set.union,
        ),
    ),
    af.AnnotDef(
        name = 'transmembrane_go',
        source = 'GO_Intercell',
        args = {
            'mainclass': 'transmembrane',
        },
    ),
    af.AnnotDef(
        name = 'transmembrane_cellphonedb',
        source = 'CellPhoneDB',
        args = {
            'transmembrane': True,
        },
    ),
    af.AnnotDef(
        name = 'transmembrane_opm',
        source = 'OPM',
        args = {
            'transmembrane': True,
        },
    ),
    af.AnnotDef(
        name = 'transmembrane_topdb',
        source = 'TopDB',
        args = {
            'topology': 'Membrane',
        },
    ),
    af.AnnotDef(
        name = 'transmembrane_locate',
        source = 'Locate',
        args = {
            'cls': {
                'typeI',
                'typeII',
                'mtmp',
            },
        },
    ),
    # adhesion
    af.AnnotDef(
        name = 'adhesion',
        source = af.AnnotOp(
            annots = (
                'adhesion_cellphonedb',
                'adhesion_go',
                'adhesion_matrisome',
                'adhesion_hgnc',
                'adhesion_integrins',
                'adhesion_zhong2015',
                'adhesion_adhesome',
            ),
            op = set.union,
        ),
    ),
    af.AnnotDef(
        name = 'adhesion_cellphonedb',
        source = 'CellPhoneDB',
        args = {
            'adhesion': bool,
        },
    ),
    af.AnnotDef(
        name = 'adhesion_integrins',
        source = 'Integrins',
    ),
    af.AnnotDef(
        name = 'adhesion_zhong2015',
        source = 'Zhong2015',
    ),
    af.AnnotDef(
        name = 'adhesion_hgnc',
        source = 'HGNC',
        args = {
            'mainclass': {
                'Type I classical cadherins',
                'Type II classical cadherins',
                '7D cadherins',
                'Desmosomal cadherins',
                'CELSR cadherins',
                'Clustered protocadherins',
                'Non-clustered protocadherins',
                'Cadherin related',
                'Integrin beta subunits',
                'Integrin alpha subunits',
                'Sialic acid binding Ig like lectins',
                'IgLON cell adhesion molecules',
                'IgCAM CXADR-related subfamily',
                'Nectins and nectin-like molecules',
                'Neurexins',
                'Neuroligins',
                (
                    'Carcinoemryonic antigen related '
                    'cell adhesion molecule family'
                ),
            }
        },
    ),
    af.AnnotDef(
        name = 'adhesion_go',
        source = 'GO_Intercell',
        args = {
            'mainclass': {'adhesion to matrix', 'adhesion to other cells'},
        }
    ),
    af.AnnotDef(
        name = 'adhesion_matrisome',
        source = af.AnnotOp(
            annots = (
                af.AnnotDef(
                    name = 'ecm_matrisome_affiliated',
                    source = 'Matrisome',
                    args = {
                        'mainclass': 'Matrisome-associated',
                        'subclass': 'ECM-affiliated Proteins',
                    },
                ),
                'cell_surface',
            ),
            op = set.intersection,
        ),
    ),
    af.AnnotDef(
        name = 'adhesion_adhesome',
        source = 'Adhesome',
    ),
    # surface enzyme
    af.AnnotDef(
        name = 'surface_enzyme',
        source = af.AnnotOp(
            annots = (
                'surface_enzyme_go',
                'surface_enzyme_surfaceome',
            ),
            op = set.union,
        ),
    ),
    af.AnnotDef(
        name = 'surface_enzyme_go',
        source = af.AnnotOp(
            annots = (
                af.AnnotOp(
                    annots = (
                        'cell_surface',
                        af.AnnotDef(
                            name = 'enzyme',
                            source = 'GO_Intercell',
                            args = {
                                'mainclass': 'enzyme',
                            },
                        ),
                        'cell_surface',
                    ),
                    op = set.intersection,
                ),
                'receptor',
            ),
            op = set.difference,
        ),
    ),
    af.AnnotDef(
        name = 'surface_enzyme_surfaceome',
        source = 'Surfaceome',
        args = {
            'mainclass': 'Enzymes',
        },
    ),
    # surface ligand
    af.AnnotDef(
        name = 'surface_ligand',
        source = af.AnnotOp(
            annots = (
                'cell_surface',
                'ligand_go',
            ),
            op = set.intersection,
        ),
    ),
    # transporter
    af.AnnotDef(
        name = 'transporter',
        source = af.AnnotOp(
            annots = (
                'transporter_surfaceome',
                'transporter_go',
            ),
            op = set.union,
        ),
    ),
    af.AnnotDef(
        name = 'transporter_surfaceome',
        source = 'Surfaceome',
        args = {
            'mainclass': 'Transporters',
        },
    ),
    af.AnnotDef(
        name = 'transporter_go',
        source = 'GO_Intercell',
        args = {
            'mainclass': {'transport', 'ion channels'},
        },
    ),
    # extracellular enzyme
    af.AnnotDef(
        name = 'extracellular_enzyme',
        source = af.AnnotOp(
            annots = (
                af.AnnotOp(
                    annots = (
                        'extracellular',
                        af.AnnotDef(
                            name = 'enzyme',
                            source = 'GO_Intercell',
                            args = {
                                'mainclass': 'enzyme',
                            },
                        ),
                    ),
                    op = set.intersection,
                ),
                af.AnnotDef(
                    name = 'matrisome_regulators',
                    source = 'Matrisome',
                    args = {
                        'subclass': 'ECM Regulators',
                    },
                ),
            ),
            op = set.union,
        ),
    ),
    # extracellular peptidase
    af.AnnotDef(
        name = 'extracellular_peptidase',
        source = af.AnnotOp(
            annots = (
                'extracellular',
                af.AnnotDef(
                    name = 'peptidase',
                    source = 'GO_Intercell',
                    args = {
                        'mainclass': 'peptidase',
                    },
                ),
            ),
            op = set.intersection,
        ),
    ),
    # growth factor binder or regulator
    af.AnnotDef(
        name = 'growth_factor_binder',
        source = 'GO_Intercell',
        args = {
            'mainclass': 'growth factor binding',
        },
    ),
    af.AnnotDef(
        name = 'growth_factor_regulator',
        source = af.AnnotOp(
            annots = (
                af.AnnotDef(
                    name = 'secreted_matrisome',
                    source = 'Matrisome',
                    args = {
                        'mainclass': 'Matrisome-associated',
                        'subclass': 'Secreted Factors',
                    },
                ),
                af.AnnotOp(
                    annots = (
                        'growth_factor_binder',
                        'extracellular_enzyme',
                    ),
                    op = set.union,
                ),
            ),
            op = set.intersection,
        ),
    ),
    # secreted
    af.AnnotDef(
        name = 'secreted',
        source = af.AnnotOp(
            annots = (
                'secreted_locate',
                'secreted_matrisome',
                'secreted_cellphonedb',
            ),
            op = set.union,
        ),
    ),
    af.AnnotDef(
        name = 'secreted_locate',
        source = 'Locate',
        args = {
            'location': 'secreted',
        }
    ),
    af.AnnotDef(
        name = 'secreted_matrisome',
        source = 'Matrisome',
        args = {
            'mainclass': 'Matrisome-associated',
            'subclass': 'Secreted Factors',
        },
    ),
    af.AnnotDef(
        name = 'secreted_cellphonedb',
        source = 'CellPhoneDB',
        args = {
            'secreted': bool,
        },
    ),
    # junctions
    af.AnnotDef(
        name = 'gap_junction',
        source = 'GO_Intercell',
        args = {
            'mainclass': 'gap junction',
        },
    ),
    af.AnnotDef(
        name = 'tight_junction',
        source = 'GO_Intercell',
        args = {
            'mainclass': 'tight junction',
        },
    ),
)

class_types = {
    'above_main': {
        'cell_surface',
        'extracellular',
        'secreted',
        'transmembrane',
        'intracellular',
    },
    'main': {
        'adhesion',
        'ecm',
        'ligand',
        'receptor',
        'surface_enzyme',
        'surface_ligand',
        'transporter',
        'extracellular_enzyme',
    },
    'small_main': {
        'gap_junction',
        'growth_factor_binder',
        'growth_factor_regulator',
        'tight_junction',
    },
    'misc': {
        'extracellular_peptidase',
        'interleukin_receptors_hgnc',
        'interleukins_hgnc',
        'chemokine_ligands_hgnc',
        'endogenous_ligands_hgnc',
    },
}


class_labels = {
    'ecm': 'Extracellular matrix',
    'interleukins_hgnc': 'Interleukins (HGNC)',
    'chemokine_ligands_hgnc': 'Chemokine ligands (HGNC)',
    'endogenous_ligands_hgnc': 'Endogenous ligands (HGNC)',
    'interleukin_receptors_hgnc': 'Interleukin receptors (HGNC)',
}


resource_labels = {
    'cellphonedb': 'CellPhoneDB',
    'topdb': 'TopDB',
    'locate': 'LOCATE',
    'go': 'Gene Ontology',
    'matrixdb': 'MatrixDB',
    'opm': 'OPM',
    'zhong2015': 'Zhong 2015',
    'kirouac': 'Kirouac 2010',
    'hgnc': 'HGNC',
    'hpmr': 'HPMR',
    'cspa': 'CSPA',
    'comppi': 'ComPPI',
    'ramilowski': 'Ramilowski 2015',
    'guide2pharma': 'Guide to Pharm',
}


def get_label(key, exceptions):
    
    return (
        exceptions[key]
            if key in exceptions else
        key.replace('_', ' ').capitalize()
    )


def get_class_label(class_key):
    
    return get_label(class_key, class_labels)


def get_resource_label(resource_key):
    
    return get_label(resource_key, resource_labels)
