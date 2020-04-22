#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#  Annotation information for intercell.py
#
#  Copyright
#  2014-2020
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#                  Nicolàs Palacio
#                  Olga Ivanova
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

import importlib as imp

import pypath.utils.go as go
import pypath.inputs.main as dataio
import pypath.internals.annot_formats as af


#TODO should go to jsons
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
        integral component of membrane OR
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
    ('LOCATE', 'extracellular'),
    ('LOCATE', 'extracellular region'),
    ('LOCATE', 'plasma membrane'),
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
        resource = af.AnnotOp(
            annots = '~receptor',
            op = set.union,
        ),
    ),
    af.AnnotDef(
        name = 'receptor_italk',
        resource = 'iTALK',
        args = {
            'mainclass': 'receptor',
        },
    ),
    af.AnnotDef(
        name = 'receptor',
        resource = 'Almen2009',
        args = {
            'mainclass': 'Receptors',
        },
    ),
    af.AnnotDef(
        name = 'receptor_cellcellinteractions',
        resource = 'CellCellInteractions',
        args = {
            'mainclass': 'Receptor',
        },
    ),
    af.AnnotDef(
        name = 'receptor_embrace',
        resource = 'EMBRACE',
        args = {
            'mainclass': 'receptor',
        },
    ),
    af.AnnotDef(
        name = 'receptor',
        resource = af.AnnotOp(
            annots = (
                'interleukin_receptor_hgnc',
            ),
            op = set.union,
        ),
    ),
    af.AnnotDef(
        name = 'interleukin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Interleukin receptors',
        },
    ),
    af.AnnotDef(
        name = 'receptor_cellphonedb',
        resource = 'CellPhoneDB',
        args = {
            'receptor': bool,
            'transmembrane': True,
        },
    ),
    af.AnnotDef(
        name = 'receptor_go',
        resource = 'GO_Intercell',
        args = {
            'mainclass': 'receptors',
        },
    ),
    af.AnnotDef(
        name = 'receptor_hpmr',
        resource = 'HPMR',
        args = {
            'role': 'Receptor',
        },
    ),
    af.AnnotDef(
        name = 'receptor_surfaceome',
        resource = 'Surfaceome',
        args = {
            'mainclass': 'Receptors',
        },
    ),
    af.AnnotDef(
        name = 'receptor_ramilowski',
        resource = 'Ramilowski2015',
        args = {
            'mainclass': 'receptor',
        },
    ),
    af.AnnotDef(
        name = 'receptor_kirouac',
        resource = 'Kirouac2010',
        args = {
            'mainclass': 'receptor',
        },
    ),
    af.AnnotDef(
        name = 'receptor_guide2pharma',
        resource = 'Guide2Pharma',
        args = {
            'mainclass': 'receptor',
        },
    ),
    af.AnnotDef(
        name = 'receptor_dgidb',
        resource = 'DGIdb',
        args = {
            'category': 'G PROTEIN COUPLED RECEPTOR',
        },
    ),
    af.AnnotDef(
        name = 'receptor_lrdb',
        resource = 'LRdb',
        args = {
            'role': 'receptor',
            'references': bool,
        },
    ),
    af.AnnotDef(
        name = 'receptor_baccin',
        resource = 'Baccin2019',
        args = {
            'mainclass': 'receptor',
        },
    ),
    af.AnnotDef(
        name = 'receptor_signalink',
        resource = 'SignaLink_function',
        args = {
            'function': 'Receptor',
        },
    ),
    # receptor subclasses from HGNC
    af.AnnotDef(
        name = 'immunoglobulin_like',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Activating leukocyte immunoglobulin like receptors',
        },
    ),
    af.AnnotDef(
        name = 'adiponectin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Adiponectin receptors',
        },
    ),
    af.AnnotDef(
        name = 'adrenalin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Adrenoceptors',
        },
    ),
    af.AnnotDef(
        name = 'angiotensin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Angiotensin receptors',
        },
    ),
    af.AnnotDef(
        name = 'vasopressin_oxytocin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Arginine vasopressin and oxytocin receptors',
        },
    ),
    af.AnnotDef(
        name = 'atypical_chemokine',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Atypical chemokine receptors',
        },
    ),
    af.AnnotDef(
        name = 'basigin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Basigin family',
        },
    ),
    af.AnnotDef(
        name = 'bombesin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Bombesin receptors',
        },
    ),
    af.AnnotDef(
        name = 'bradykinin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Bradykinin receptors',
        },
    ),
    af.AnnotDef(
        name = 'butyrophilin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Butyrophilins',
        },
    ),
    af.AnnotDef(
        name = 'cc_chemokine',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'C-C motif chemokine receptors',
        },
    ),
    af.AnnotDef(
        name = 'cx3c_chemokine',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'C-X-3-C motif chemokine receptors',
        },
    ),
    af.AnnotDef(
        name = 'cxc_chemokine',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'C-X-C motif chemokine receptors',
        },
    ),
    af.AnnotDef(
        name = 'celsr_cadherin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'CELSR cadherins',
        },
    ),
    af.AnnotDef(
        name = '5ht_gprotein',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': '5-hydroxytryptamine receptors, G protein-coupled',
        },
    ),
    af.AnnotDef(
        name = '5ht_ionotropic',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': '5-hydroxytryptamine receptors, ionotropic',
        },
    ),
    af.AnnotDef(
        name = 'activating_leukocyte_ig_like',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Activating leukocyte immunoglobulin like receptors',
        },
    ),
    af.AnnotDef(
        name = 'adenosine',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Adenosine receptors',
        },
    ),
    af.AnnotDef(
        name = 'calcitonin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Calcitonin receptors',
        },
    ),
    af.AnnotDef(
        name = 'calcium',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Calcium sensing receptors',
        },
    ),
    af.AnnotDef(
        name = 'cannabinoid',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Cannabinoid receptors',
        },
    ),
    af.AnnotDef(
        name = 'chemerin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Chemerin receptor',
        },
    ),
    af.AnnotDef(
        name = 'cholecystokinin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Cholecystokinin receptors',
        },
    ),
    af.AnnotDef(
        name = 'muscarinic_cholinergic',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Cholinergic receptors muscarinic',
        },
    ),
    af.AnnotDef(
        name = 'nicotinic_cholinergic',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Cholinergic receptors nicotinic subunits',
        },
    ),
    af.AnnotDef(
        name = 'collectin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Collectins',
        },
    ), # innate immunity receptors for sugar and lipid patterns
    af.AnnotDef(
        name = 'complement_gpcr',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Complement component GPCRs',
        },
    ), # receptors for chemotactic immune signals
    af.AnnotDef(
        name = 'crh',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Corticotropin releasing hormone receptors',
        },
    ),
    af.AnnotDef(
        name = 'dopamine',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Dopamine receptors',
        },
    ),
    af.AnnotDef(
        name = 'ephrin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'EPH receptors',
        },
    ),
    af.AnnotDef(
        name = 'endothelin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Endothelin receptors',
        },
    ),
    af.AnnotDef(
        name = 'erbb_rtk',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Erb-b2 receptor tyrosine kinases',
        },
    ),
    af.AnnotDef(
        name = 'f2r',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'F2R receptors',
        },
    ), # GPCRs for thrombin and trypsin
    af.AnnotDef(
        name = 'formyl_peptide',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Formyl peptide receptors',
        },
    ), # formyl-methionyl peptides are neutrophil chemoattractants
    af.AnnotDef(
        name = 'free_fatty_acid',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Free fatty acid receptors',
        },
    ),  # intestinal short chain fatty acid GPCRs, regulating
        # whole-body energy homeostasis
    af.AnnotDef(
        name = 'bile_acid',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'G protein-coupled bile acid receptor',
        },
    ),  # GPCR for bile acid
    af.AnnotDef(
        name = 'estrogen_gpcr',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'G protein-coupled estrogen receptor',
        },
    ),  # although this receptor is intracellular, there is no reason we
        # shouldn't treat it the same way as plasma membrane receptors
    af.AnnotDef(
        name = 'nuclear_hormone',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Nuclear hormone receptors',
        },
    ),  # although these receptors are intracellular, there is no reason we
        # shouldn't treat it the same way as plasma membrane receptors
    af.AnnotDef(
        name = 'gpcr_orphan',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': {
                'G protein-coupled receptors, Class A orphans',
                'G protein-coupled receptors, Class C orphans',
            },
        },
    ),  # GPCRs mostly without known ligand, all in the cell membrane
    af.AnnotDef(
        name = 'frizzled',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'G protein-coupled receptors, Class F frizzled',
        },
    ),  # GPCRs for Wnt
    af.AnnotDef(
        name = 'galanin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Galanin receptors',
        },
    ),
    af.AnnotDef(
        name = 'gaba',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': {
                'Gamma-aminobutyric acid type A receptor subunits',
                'Gamma-aminobutyric acid type B receptor subunits',
            }
        },
    ),
    af.AnnotDef(
        name = 'glucagon',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Glucagon receptor family',
        },
    ),
    af.AnnotDef(
        name = 'glutamate_ionotropic',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': {
                'Glutamate ionotropic receptor AMPA type subunits',
                'Glutamate ionotropic receptor NMDA type subunits',
                'Glutamate ionotropic receptor delta type subunits',
                'Glutamate ionotropic receptor kainate type subunits',
            }
        },
    ),
    af.AnnotDef(
        name = 'glutamate_metabotropic',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Glutamate metabotropic receptors',
        },
    ),
    af.AnnotDef(
        name = 'glycine',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Glycine receptors',
        },
    ),
    af.AnnotDef(
        name = 'glycoprotein_hormone',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Glycoprotein hormone receptors',
        },
    ),  # receptors for TSH, FSH and LH
    af.AnnotDef(
        name = 'gonadotropin_releasing_hormone',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Gonadotropin releasing hormone receptors',
        },
    ),  # receptors for GnRH
    af.AnnotDef(
        name = 'natriuretic_peptide',
        parent = 'receptor',
        resource = {'P16066', 'P20594'},
    ),
    af.AnnotDef(
        name = 'guanilyn',
        parent = 'receptor',
        resource = {'P25092'},
    ),
    af.AnnotDef(
        name = 'histamine',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Histamine receptors',
        },
    ),
    af.AnnotDef(
        name = 'hydroxycarboxylic_acid',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Hydroxy-carboxylic acid receptors',
        },
    ),  # receptors for lactate, niacin, etc
    af.AnnotDef(
        name = 'orexin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Hypocretin receptors',
        },
    ),
    af.AnnotDef(
        name = 'inhibitory_leukocyte_ig_like',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Inhibitory leukocyte immunoglobulin like receptors',
        },
    ),
    af.AnnotDef(
        name = 'interferon',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Interferon receptors',
        },
    ),
    af.AnnotDef(
        name = 'interleukin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Interleukin receptors',
        },
    ),
    af.AnnotDef(
        name = 'killer_cell_ig_like',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Killer cell immunoglobulin like receptors',
        },
    ),  # receptors for HLAs
    af.AnnotDef(
        name = 'killer_cell_lectin_like',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Killer cell lectin like receptors',
        },
    ),  # receptors mostly for HLAs
    af.AnnotDef(
        name = 'ly6_plaur',
        parent = 'receptor',
        resource = {
            'Q03405', 'Q8IV16',
        }
    ),
    af.AnnotDef(
        name = 'ly6_plaur',
        parent = 'receptor_regulator',
        resource = {
            'Q5SQ64', 'Q8N2G4', 'Q86Y78', 'Q8N6Q3', 'Q16553', 'P0DP58',
            'O94772', 'P13987', 'P0DP57', 'O43653', 'Q8NI32', 'P0C8F1',
        },
    ),  #
    af.AnnotDef(
        name = 'leukotriene',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Leukotriene receptors',
        },
    ),  # receptors for leukotrienes, eicosanoids, N-formyl-met peptides, etc
    af.AnnotDef(
        name = 'pentraxin',
        parent = 'receptor',
        resource = {'Q15818', 'O95502'},
    ),
    af.AnnotDef(
        name = 'short_pentraxin',
        parent = 'secreted_receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Short pentraxins',
        },
    ),  # bind to microbial antigens, DNA and histones
    af.AnnotDef(
        name = 'low_density_lipoprotein',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Low density lipoprotein receptors',
        },
    ),  # receptors for low density lipoproteins
    af.AnnotDef(
        name = 'lysophosphatidic_acid',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Lysophosphatidic acid receptors',
        },
    ),
    af.AnnotDef(
        name = 'melanin_concentrating_hormone',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Melanin concentrating hormone receptors',
        },
    ),
    af.AnnotDef(
        name = 'melanocortin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Melanocortin receptors',
        },
    ),
    af.AnnotDef(
        name = 'melatonin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Melatonin receptors',
        },
    ),
    af.AnnotDef(
        name = 'membrane_associated_progesterone',
        parent = 'receptor',
        resource = {'O15173', 'O00264'},
    ),  # PGRMC1 is in SER and microsome membrane, but it does not matter
    af.AnnotDef(
        name = 'neuromedin_u',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Neuromedin U receptors',
        },
    ),
    af.AnnotDef(
        name = 'neuropeptide',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Neuropeptide receptors',
        },
    ),
    af.AnnotDef(
        name = 'neurotensin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Neurotensin receptors',
        },
    ),
    af.AnnotDef(
        name = 'notch',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Notch receptors',
        },
    ),
    af.AnnotDef(
        name = 'olfactory',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': {
                'Olfactory receptors, family 1',
                'Olfactory receptors, family 2',
                'Olfactory receptors, family 3',
                'Olfactory receptors, family 4',
                'Olfactory receptors, family 5',
                'Olfactory receptors, family 6',
                'Olfactory receptors, family 7',
                'Olfactory receptors, family 8',
                'Olfactory receptors, family 9',
                'Olfactory receptors, family 10',
                'Olfactory receptors, family 11',
                'Olfactory receptors, family 12',
                'Olfactory receptors, family 13',
                'Olfactory receptors, family 14',
                'Olfactory receptors, family 51',
                'Olfactory receptors, family 52',
                'Olfactory receptors, family 56',
            },
        },
    ),  # odorant receptors
    af.AnnotDef(
        name = 'opioid',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Opioid receptors',
        },
    ),
    af.AnnotDef(
        name = 'opsin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Opsin receptors',
        },
    ),
    af.AnnotDef(
        name = 'oxoglutarate',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Oxoglutarate receptor',
        },
    ),
    af.AnnotDef(
        name = 'p2y_purinergic',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'P2Y receptors',
        },
    ),  # receptors for ADP, ATP, UDP, UTP
    af.AnnotDef(
        name = 'p2x_purinergic',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Purinergic receptors P2X',
        },
    ),  # receptors for ATP
    af.AnnotDef(
        name = 'parathyroid_hormone',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Parathyroid hormone receptors',
        },
    ),
    af.AnnotDef(
        name = 'peptide',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Peptide receptors',
        },
    ),  # receptors for TRH, motilin, apelin, PrRP, QRFP, ghrelin, etc
    af.AnnotDef(
        name = 'platelet_activating_factor',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Platelet activating factor receptor',
        },
    ),
    af.AnnotDef(
        name = 'plexin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Plexins',
        },
    ),  # receptors for semaphorins
    af.AnnotDef(
        name = 'progestin_and_adipoq',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Progestin and adipoQ receptor family',
        },
        exclude = {'Q6TCH7', 'Q8IY49', 'Q15546'},
    ),  # G protein coupled progesterone and ADIPOQ hormone receptors
    af.AnnotDef(
        name = 'prokineticin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Prokineticin receptors',
        },
    ),
    af.AnnotDef(
        name = 'prostaglandin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Prostaglandin receptors',
        },
    ),
    af.AnnotDef(
        name = 'receptor_tyrosine_phosphatase',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Protein tyrosine phosphatases receptor type',
        },
        exclude = {'Q92932'},
    ),
    af.AnnotDef(
        name = 'proteoglycan',
        parent = 'receptor',
        resource = {'P16070', 'Q6UVK1'},
    ),
    af.AnnotDef(
        name = 'proteoglycan',
        parent = 'receptor_regulator',
        resource = {'O00468'},
    ),
    af.AnnotDef(
        name = 'pseudoautosomal_region',
        parent = 'receptor',
        resource = {'P15509', 'Q86VZ1', 'Q9HC73', 'P26951', 'Q01113'},
    ),
    af.AnnotDef(
        name = 'relt',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'RELT family',
        },
    ),
    af.AnnotDef(
        name = 'gpcr_activity_modifying',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': (
                'Receptor (G protein-coupled) activity modifying protein'
            ),
        },
    ),  # receptors for adrenomedullin
    af.AnnotDef(
        name = 'receptor_transporter',
        parent = 'receptor_regulator',
        resource = 'HGNC',
        args = {
            'mainclass': 'Receptor transporter proteins',
        },
    ),  # regulate GPCRs, especially taste and olfactory receptors
    af.AnnotDef(
        name = 'tyrosine_kinase',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Receptor tyrosine kinases',
        },
    ),
    af.AnnotDef(
        name = 'relaxin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Relaxin family peptide receptors',
        },
    ),
    af.AnnotDef(
        name = 'repulsive_guidance',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Repulsive guidance molecule family',
        },
    ),  # BMP coreceptors
    af.AnnotDef(
        name = 'slitrk',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'SLIT and NTRK like family',
        },
    ),  # receptors regulating synapse development in CNS
    af.AnnotDef(
        name = 'scavenger_receptor_cysteine_rich',
        parent = 'receptor',
        resource = {'Q9UEW3', 'Q6ZMJ2', 'P21757', 'P06127', 'Q86VB7'},
    ),
    af.AnnotDef(
        name = 'scavenger_receptor_cysteine_rich',
        parent = 'secreted_receptor',
        resource = {'Q9UGM3', 'A1L4H1', 'Q86VB7'},
    ),

    # secreted receptors
    af.AnnotDef(
        name = 'ly6_plur',
        parent = 'secreted_receptor',
        resource = {'Q6UX82', 'P55000', 'P13987'},
    ),
    af.AnnotDef(
        name = 'pentraxin',
        parent = 'secreted_receptor',
        resource = {'Q96A99', 'P26022', 'P47972'},
    ),
    af.AnnotDef(
        name = 'ms4',
        parent = 'receptor_regulator',
        resource = af.AnnotOp(
            annots = (
                af.AnnotDef(
                    name = 'ms4',
                    resource = 'HGNC',
                    args = {
                        'mainclass': 'Membrane spanning 4-domains',
                    },
                ),
                'plasma_membrane',
            ),
            op = set.intersection,
        ),
    ),
    af.AnnotDef(
        name = 'ige',
        parent = 'receptor',
        resource = {'Q01362'},
    ),
    af.AnnotDef(
        name = 'peptidoglycan',
        parent = 'secreted_receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Peptidoglycan recognition proteins',
        },
    ),  # apart from peptide recognition they have anti-microbial activity
        # either enzymatically or other ways
    af.AnnotDef(
        name = 'scavenger',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Scavenger receptors',
        },
        exclude = {'Q8WTU2', 'Q6AZY7', 'A1L4H1', 'Q14108'}
    ),
    af.AnnotDef(
        name = 'scavenger',
        parent = 'secreted_receptor',
        resource = {'Q8WTU2', 'A1L4H1', 'Q86VB7'},
    ),
    af.AnnotDef(
        name = 'sialic_acid_binding_lectin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Sialic acid binding Ig like lectins',
        },
    ),  # not all of them are receptors, some of them might be cell-cell
        # adhesion proteins, depending on the intracellular domains
    af.AnnotDef(
        name = 'somatostatin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Somatostatin receptors',
        },
    ),
    af.AnnotDef(
        name = 'sphingosine_phosphate',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Sphingosine 1-phosphate receptors',
        },
    ),
    af.AnnotDef(
        name = 'succinate',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Succinate receptor',
        },
    ),
    af.AnnotDef(
        name = 't_cell',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': {
                'T cell receptor alpha locus at 14q11.2',
                'T cell receptor beta locus at 7q34',
                'T cell receptor delta locus at 14q11.2',
                'T cell receptor gamma locus at 7p14',
            },
        },
    ),
    af.AnnotDef(
        name = 'tir_domain',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'TIR domain containing',
        },
        exclude = {'Q8IUC6', 'Q99836', 'Q6SZW1', 'Q86XR7'}
    ),
    af.AnnotDef(
        name = 'tachykinin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Tachykinin receptors',
        },
    ),
    af.AnnotDef(
        name = 'taste',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': {
                'Taste 1 receptors',
                'Taste 2 receptors',
            },
        },
    ),
    af.AnnotDef(
        name = 'toll_like',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Toll like receptors',
        },
    ),
    af.AnnotDef(
        name = 'trace_amin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Trace amine receptors',
        },
    ),
    af.AnnotDef(
        name = 'tumor_necrosis_factor',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Tumor necrosis factor receptor superfamily',
        },
        exclude = {'O9540', 'O0030'},
    ),
    af.AnnotDef(
        name = 'tumor_necrosis_factor',
        parent = 'secreted_receptor',
        resource = {'O9540', 'O0030'},
    ),
    af.AnnotDef(
        name = 'type1_serine_threonine_kinase',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Type 1 receptor serine/threonine kinases',
        },
    ),
    af.AnnotDef(
        name = 'type2_serine_threonine_kinase',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Type 2 receptor serine/threonine kinases',
        },
    ),
    af.AnnotDef(
        name = 'vasoactive_intestinal_peptide',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Vasoactive intestinal peptide receptor family',
        },
    ),
    af.AnnotDef(
        name = 'vomeronasal',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Vomeronasal receptors',
        },
    ),
    af.AnnotDef(
        name = 'xc_motif_chemokine',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'X-C motif chemokine receptors',
        },
    ),
    # subclasses from Almen 2009
    af.AnnotDef(
        name = 'transforming_growth_factor',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'Act.TGFB',
        },
    ),
    af.AnnotDef(
        name = 'axl',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'Axl',
        },
    ),  # their ligands are in the ECM
    af.AnnotDef(
        name = 'butyrophilin',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'Butyrophilin',
        },
    ),  # not clear if these are all receptors or some of them are ligands,
        # transporters or regulators of other membrane proteins
    af.AnnotDef(
        name = 'scavenger',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'SCAR',
        },
    ),
    af.AnnotDef(
        name = 'cytokine',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'CytokineR',
        },
    ),
    af.AnnotDef(
        name = 'egf',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'EGFR',
        },
    ),
    af.AnnotDef(
        name = 'ephrin',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'Eph',
        },
    ),
    af.AnnotDef(
        name = 'fgf',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'FGFR',
        },
    ),
    af.AnnotDef(
        name = 'fc',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'FcR',
        },
    ),
    af.AnnotDef(
        name = 'frizzled',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'Frizzled',
        },
    ),
    af.AnnotDef(
        name = 'gpcr',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'GPCR',
        },
    ),
    af.AnnotDef(
        name = 'glutamate',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'Glutamate',
        },
    ),
    af.AnnotDef(
        name = 'ig_like',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'mainclass': 'Receptor',
            'classes': 'IG',
        },
    ),
    af.AnnotDef(
        name = 'il17',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'mainclass': 'Receptor',
            'classes': 'IL17',
        },
    ),
    af.AnnotDef(
        name = 'igf',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'InsR',
        },
    ),
    af.AnnotDef(
        name = 'killer_cell_ig_like',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'KIR',
        },
    ),
    af.AnnotDef(
        name = 'kinase',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'mainclass': 'Receptors',
            'classes': {
                'Kinase',
                'KInase' # this is a typo
            },
        },
    ),  # all receptors with kinase activity,
        # I think it's useful to have such category
    af.AnnotDef(
        name = 'low_density_lipoprotein',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'LDLR',
        },
    ),
    af.AnnotDef(
        name = 'leukocyte_ig_like',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'LILR',
        },
    ),
    af.AnnotDef(
        name = 'mannose',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'MacrophageMannoseR',
        },
    ),
    af.AnnotDef(
        name = 'netrin',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'NetrinR',
        },
    ),
    af.AnnotDef(
        name = 'neuropilin',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'Neuropilin',
        },
    ),  # receptor for semaphorins, VEGF and PLGF
    af.AnnotDef(
        name = 'notch',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'Notch',
        },
    ),
    af.AnnotDef(
        name = 'olfactory',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': {'Olf', 'PutativeOlfR'},
        },
    ),
    af.AnnotDef(
        name = 'cd1',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'OtherCD1',
        },
    ),
    af.AnnotDef(
        name = 'cd300',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'OtherCD300',
        },
    ),
    af.AnnotDef(
        name = 'natural_cytotoxicity_triggering',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'OtherNCR',
        },
    ),
    af.AnnotDef(
        name = 'poliovirus',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'OtherPVR',
        },
    ),
    af.AnnotDef(
        name = 'roundabout',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'OtherROBO',
        },
    ),
    af.AnnotDef(
        name = 'triggering',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'OtherTREM',
        },
        exclude = {'Q6UXN2'}
    ),
    af.AnnotDef(
        name = 'pdgf',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'PDGFR',
        },
    ),
    af.AnnotDef(
        name = 'patched',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'Patched',
        },
    ),
    af.AnnotDef(
        name = 'plexin',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'Plexin',
        },
    ),
    af.AnnotDef(
        name = 'receptor_activity_modifying',
        parent = 'receptor_regulator',
        resource = 'Almen2009',
        args = {
            'classes': 'RAMP',
        },
    ),  # affect other receptors such as CALCRL
    af.AnnotDef(
        name = 'retinal_guanylyl_cyclase',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'RGC',
        },
    ),  # receptors for various compounds, e.g. natriuretic peptide
        # and E.coli enterotoxin
    af.AnnotDef(
        name = 'receptor_transporter',
        parent = 'receptor_regulator',
        resource = 'Almen2009',
        args = {
            'classes': 'RTP',
        },
    ),  # regulate other receptors,
        # especially their trafficking to the plasma membrane
    af.AnnotDef(
        name = 'receptor_tyrosine_phosphatase',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'ReceptorTypePhosphatases',
        },
    ),  # not sure all these are receptors
    af.AnnotDef(
        name = 'rhodopsin',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'Rhodopsin',
        },
    ),  # GPCRs for various ligands, hormones, neurotransmitters, etc
    af.AnnotDef(
        name = 'secretin',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'Secretin',
        },
    ),  # GPCRs for various hormones
    af.AnnotDef(
        name = 'syndecan',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'Syndecan',
        },
    ),  # heparan sulphate carrying cell surface proteins
        # transferring signals to the cytoskeleton, regulating cell shape
    af.AnnotDef(
        name = 'taste',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'TAS2R',
        },
    ),
    af.AnnotDef(
        name = 't_cell',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'TCR',
        },
    ),
    af.AnnotDef(
        name = 'tumor_necrosis_factor',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'TNFNGF',
        },
    ),  # maybe this is a broader superfamily
    af.AnnotDef(
        name = 'toll_like',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'TOLL',
        },
    ),  # receptors for various pathogen patterns e.g. LPS
    af.AnnotDef(
        name = 'teneurin',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'Teneurin',
        },
    ),  # functioning in neural development and maybe elsewhere;
        # these are definitely receptors, but also ligands and
        # maybe cell-cell adhesion molecules
    af.AnnotDef(
        name = 'transferrin',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'Transferrin',
        },
    ),
    af.AnnotDef(
        name = 'type1_ig_like_cytokine',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'Type1',
        },
    ),  # mostly interleukin receptors
    af.AnnotDef(
        name = 'type2_ig_like_cytokine',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'Type2',
        },
    ),  # mostly interferon and interleukin receptors
    af.AnnotDef(
        name = 'vomeronasal',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'V1R',
        },
    ),  # pheromone receptors
    af.AnnotDef(
        name = 'gaba_ach',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'cys-loop',
        },
    ),  # receptors for ACh and GABA
    af.AnnotDef(
        name = 'neurotrophin',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'neutrophin',
        },
    ),

    # ECM
    af.AnnotDef(
        name = 'ecm',
        resource = af.AnnotOp(
            annots = (
                'ecm_matrixdb',
                'ecm_matrisome',
                'ecm_go',
                'ecm_ramilowski',
                'ecm_cellcellinteractions',
                'ecm_uniprot',
            ),
            op = set.union,
        ),
    ),
    af.AnnotDef(
        name = 'ecm_matrisome',
        resource = af.AnnotOp(
            annots = (
                af.AnnotDef(
                    name = 'ecm_matrisome_core',
                    resource = 'Matrisome',
                    args = {
                        'mainclass': 'Core matrisome',
                    },
                ),
                af.AnnotOp(
                    annots = (
                        af.AnnotDef(
                            name = 'ecm_matrisome_affiliated',
                            resource = 'Matrisome',
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
        resource = 'MatrixDB',
        args = {
            'mainclass': 'ecm',
        },
    ),
    af.AnnotDef(
        name = 'ecm_go',
        resource = 'GO_Intercell',
        args = {
            'mainclass': 'ecm structure',
        },
    ),
    af.AnnotDef(
        name = 'ecm_ramilowski',
        resource = 'Ramilowski_location',
        args = {
            'location': {
                'extracellular matrix',
                'basement membrane',
            },
        },
    ),
    af.AnnotDef(
        name = 'ecm_uniprot',
        resource = 'UniProt_location',
        args = {
            'location': 'Extracellular matrix',
        },
    ),
    af.AnnotDef(
        name = 'ecm_cellcellinteractions',
        resource = 'CellCellInteractions',
        args = {
            'mainclass': 'ECM',
        },
    ),
    # specific subclasses from HGNC
    af.AnnotDef(
        name = 'collagen_proteoglycan_ecm',
        resource = 'HGNC',
        args = {
            'mainclass': 'Collagen proteoglycans',
        },
    ),
    af.AnnotDef(
        name = 'collagen_ecm',
        resource = 'HGNC',
        args = {
            'mainclass': 'Collagens',
        },
    ),
    af.AnnotDef(
        name = 'emi_ecm',
        resource = 'HGNC',
        args = {
            'mainclass': 'EMI domain containing',
        },
    ),  # this could be also cell-matrix adhesion, although these proteins
        # are not in the cell membrane but all secreted
    af.AnnotDef(
        name = 'fibrillin_ecm',
        resource = 'HGNC',
        args = {
            'mainclass': 'Fibrillins',
        },
    ),
    af.AnnotDef(
        name = 'laminin_ecm',
        resource = 'HGNC',
        args = {
            'mainclass': 'Laminin subunits',
        },
    ),
    af.AnnotDef(
        name = 'fibulin_ecm',
        resource = 'HGNC',
        args = {
            'mainclass': 'Fibulins',
        },
    ),  # parts of ECM, especially elastic fibers, one of them is a ligand
        # for EGFR (but still an EVM protein at the same time)
    af.AnnotDef(
        name = 'hyalectan_proteoglycan_ecm',
        resource = 'HGNC',
        args = {
            'mainclass': 'Hyalectan proteoglycans',
        },
    ),
    af.AnnotDef(
        name = 'matrilin_ecm',
        resource = 'HGNC',
        args = {
            'mainclass': 'Matrilins',
        },
    ),  # cartilage ECM
    af.AnnotDef(
        name = 'mucin_ecm',
        resource = af.AnnotOp(
            annots = (
                af.AnnotDef(
                    name = 'mucin_hgnc',
                    source ='HGNC',
                    args = {
                        'mainclass': 'Mucins',
                    },
                ),
                'secreted'
            ),
            op = set.intersection
        ),
    ),
    af.AnnotDef(
        name = 'proteoglycan_ecm',
        resource = {'O00468', 'P98160'},
    ),
    af.AnnotDef(
        name = 'sibling_ecm',
        resource = 'HGNC',
        args = {
            'mainclass': 'SIBLING family',
        },
    ),
    af.AnnotDef(
        name = 'sparc_ecm_regulator',
        resource = 'HGNC',
        args = {
            'mainclass': 'SPARC family',
        },
    ),  # act either on ligands or ECM or both
    af.AnnotDef(
        name = 'small_leucine_rich_repeat_proteoglycan_ecm',
        resource = 'HGNC',
        args = {
            'mainclass': 'Small leucine rich repeat proteoglycans',
        },
    ),
    af.AnnotDef(
        name = 'zona_pellucida_glycoprotein_ecm',
        resource = 'HGNC',
        args = {
            'mainclass': 'Zona pellucida glycoproteins',
        },
    ),  # ECM of the zona pellucida (zone surrounding the oocyte)

    # ligand
    af.AnnotDef(
        name = 'ligand',
        resource = af.AnnotOp(
            annots = (
                'ligand_cellphonedb',
                'ligand_go',
                'ligand_hpmr',
                'ligand_ramilowski',
                'ligand_kirouac',
                'ligand_guide2pharma',
                'ligand_hgnc',
                'ligand_dgidb',
                'ligand_baccin',
                'ligand_signalink',
                'ligand_italk',
                'ligand_cellcellinteractions',
                'ligand_embrace',
            ),
            op = set.union,
        ),
    ),
    af.AnnotDef(
        name = 'ligand_italk',
        resource = 'iTALK',
        args = {
            'mainclass': 'ligand',
        },
    ),
    af.AnnotDef(
        name = 'ligand_cellcellinteractions',
        resource = 'CellCellInteractions',
        args = {
            'mainclass': 'Ligand',
        },
    ),
    af.AnnotDef(
        name = 'ligand_embrace',
        resource = 'EMBRACE',
        args = {
            'mainclass': 'ligand',
        },
    ),
    af.AnnotDef(
        name = 'ligand',
        resource = af.AnnotOp(
            annots = (
                'interleukin_hgnc',
                'endogenous_ligand_hgnc',
                'chemokine_ligand_hgnc',
            ),
            op = set.union,
        ),
    ),
    af.AnnotDef(
        name = 'ligand_cellphonedb',
        resource = 'CellPhoneDB',
        args = {
            'secreted': bool,
        },
    ),
    af.AnnotDef(
        name = 'ligand_go',
        resource = 'GO_Intercell',
        args = {
            'mainclass': 'ligands',
        },
    ),
    af.AnnotDef(
        name = 'ligand_hpmr',
        resource = 'HPMR',
        args = {
            'role': 'Ligand',
        },
    ),
    af.AnnotDef(
        name = 'ligand_ramilowski',
        resource = 'Ramilowski2015',
        args = {
            'mainclass': 'ligand',
        },
    ),
    af.AnnotDef(
        name = 'ligand_kirouac',
        resource = 'Kirouac2010',
        args = {
            'mainclass': 'ligand',
        },
    ),
    af.AnnotDef(
        name = 'ligand_guide2pharma',
        resource = 'Guide2Pharma',
        args = {
            'mainclass': 'ligand',
        },
    ),
    af.AnnotDef(
        name = 'ligand_dgidb',
        resource = af.AnnotOp(
            annots = (
                'growth_factor_dgidb',
                'hormone_dgidb',
            ),
            op = set.union,
        ),
    ),
    af.AnnotDef(
        name = 'growth_factor_dgidb',
        resource = 'DGIdb',
        args = {
            'category': 'GROWTH FACTOR',
        },
    ),
    af.AnnotDef(
        'hormone_dgidb',
        resource = 'DGIdb',
        args = {
            'category': 'HORMONE ACTIVITY',
        },
    ),
    af.AnnotDef(
        name = 'ligand_lrdb',
        resource = 'LRdb',
        args = {
            'role': 'ligand',
            'references': bool,
        },
    ),
    af.AnnotDef(
        name = 'ligand_baccin',
        resource = 'Baccin2019',
        args = {
            'mainclass': 'ligand',
        },
    ),
    af.AnnotDef(
        name = 'ligand_signalink',
        resource = 'SignaLink_function',
        args = {
            'function': 'Ligand',
        },
    ),
    # ligands from HGNC
    af.AnnotDef(
        name = 'angiopoietin',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Angiopoietin like family',
        },
    ),
    af.AnnotDef(
        name = 'basigin',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Basigin family',
        },
    ),
    af.AnnotDef(
        name = 'bmp',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Bone morphogenetic proteins',
        },
    ),
    af.AnnotDef(
        name = 'c1q',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'C1q and TNF related',
        },
    ),
    af.AnnotDef(
        name = 'ccn',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Cellular communication network factors',
        },
    ),
    af.AnnotDef(
        name = 'interleukin',
        resource = 'HGNC',
        args = {
            'mainclass': 'Interleukins',
        },
    ),
    af.AnnotDef(
        name = 'endogenous',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Endogenous ligands',
        },
    ), # a very few among these are actually not secreted
    af.AnnotDef(
        name = 'chemokine',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Chemokine ligands',
        },
    ),
    af.AnnotDef(
        name = 'neurotrophin',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Neurotrophins',
        },
    ),
    af.AnnotDef(
        name = 'gdnf',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'GDNF family ligands',
        },
    ),
    af.AnnotDef(
        name = 'chordin',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Chordin family',
        },
    ), # BMP antagonist ligands
    af.AnnotDef(
        name = 'cysteine_rich_bmp_regulator',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Cysteine rich transmembrane BMP regulators',
        },
    ), # BMP agonist and antagonist ligands
    af.AnnotDef(
        name = 'dan',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'DAN family',
        },
    ), # TGF & BMP signaling agonists and antagonists
    af.AnnotDef(
        name = 'fgf',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Fibroblast growth factor family',
        },
    ), # with the exception of FGF13: that's not secreted
    af.AnnotDef(
        name = 'gdnf',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'GDNF family ligands',
        },
    ), # neurotrophic ligands
    af.AnnotDef(
        name = 'growth_hormone',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Growth hormone family',
        },
    ),
    af.AnnotDef(
        name = 'hedgehog',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Hedgehog signaling molecule family',
        },
    ),  # hedgehog proteins are initially membrane bound and can be
        # solubilized later
    af.AnnotDef(
        name = 'hdgf',
        parent = 'ligand',
        resource = {'P51858'},
    ),  # hepatoma-derived growth factor (HDGFL1 is not included because
        # little is known about its role and whether it's secreted)
    af.AnnotDef(
        name = 'igf',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'IGF like family',
        },
    ),
    af.AnnotDef(
        name = 'izumo',
        parent = 'ligand',
        resource = {'Q1ZYL8'},
    ),  # ligands in sperm-egg fusion
    af.AnnotDef(
        name = 'inhibin',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Inhibin subunits',
        },
    ),
    af.AnnotDef(
        name = 'interleukin6',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Interleukin 6 type cytokine family',
        },
    ),
    af.AnnotDef(
        name = 'interleukin',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Interleukins',
        },
    ),
    af.AnnotDef(
        name = 'leucine_rich_glioma_inactivated',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'LGI family',
        },
    ),  # maybe not ligands in a strict sense but don't fit either
        # in other categories
    af.AnnotDef(
        name = 'tgf_beta_binding',
        parent = 'ligand_regulator',
        resource = 'HGNC',
        args = {
            'mainclass': (
                'Latent transforming growth factor beta binding proteins'
            ),
        },
    ),
    af.AnnotDef(
        name = 'mia',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'MIA family',
        },
        exclude = {'Q5JRA6', 'Q96PC5'},
    ),
    af.AnnotDef(
        name = 'neuferricin_neudensin',
        parent = 'ligand',
        resource = {'Q8WUJ1', 'Q9UMX5'},
    ),
    af.AnnotDef(
        name = 'netrin',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Netrins',
        },
    ),  # secreted axon guidance molecules
    af.AnnotDef(
        name = 'neurotrophin',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Neurotrophins',
        },
    ),
    af.AnnotDef(
        name = 'oocyte_secreted',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'OOSP family',
        },
    ),  # not sure these are ligands, but this category looks the most likely
        # at least in a recent paper PLAC1 has been described to activate
        # FGFR2 together with FGF7
    af.AnnotDef(
        name = 'prostate_and_testis_expressed',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'PATE family',
        },
    ),  # ligands modulating nicotinic ACh receptors and sperm motility
    af.AnnotDef(
        name = 'pregnancy_specific_glycoprotein',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Pregnancy specific glycoproteins',
        },
    ),  # bind to cell surface moieties, regulate other ligands --
        # overall they fit the best to the ligand category
    af.AnnotDef(
        name = 'proteoglycan',
        parent = 'ligand_regulator',
        resource = {'Q03167'},
    ),
    af.AnnotDef(
        name = 's100_calcium_binding',
        parent = 'ligand',
        resource = {'P31151', 'P80511', 'P05109'},
    ),
    af.AnnotDef(
        name = 'sparc',
        parent = 'ligand_regulator',
        resource = 'HGNC',
        args = {
            'mainclass': 'SPARC family',
        },
    ),
    af.AnnotDef(
        name = 'scavenger_receptor_cysteine_rich',
        parent = 'ligand_regulator',
        resource = {'Q8WTU2', ''},
    ),
    af.AnnotDef(
        name = 'scavenger_receptor_cysteine_rich',
        parent = 'ligand',
        resource = {'O43866'},
    ),
    af.AnnotDef(
        name = 'frizzled_related',
        parent = 'ligand_regulator',
        resource = 'HGNC',
        args = {
            'mainclass': 'Secreted frizzled-related proteins',
        },
    ),  # secreted proteins binding WNT ligands
    af.AnnotDef(
        name = 'secretoglobin',
        parent = 'ligand_regulator',
        resource = 'HGNC',
        args = {
            'mainclass': 'Secretoglobins',
        },
    ),  # secreted proteins binding small molecule ligands
    af.AnnotDef(
        name = 'tafa_chemokine_like',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'TAFA chemokine like family',
        },
    ),
    af.AnnotDef(
        name = 'prosalusin',
        parent = 'ligand',
        resource = {'Q8N2E6'},
    ),
    af.AnnotDef(
        name = 'transforming_growth_factor_beta',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Transforming growth factor beta family',
        },
    ),
    af.AnnotDef(
        name = 'tumor_necrosis_factor',
        parent = 'surface_ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Tumor necrosis factor superfamily',
        },
        exclude = {'O75888'},
    ),
    af.AnnotDef(
        name = 'tumor_necrosis_factor',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Tumor necrosis factor superfamily',
        },
        exclude = {
            'Q9UNG2', 'P32970', 'O14788', 'P48023',
            'P32971', 'P41273', 'P23510', 'Q06643',
        },
    ),
    af.AnnotDef(
        name = 'vascular_endothelial_growth_factor',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'VEGF family',
        },
    ),
    af.AnnotDef(
        name = 'wnt',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Wnt family',
        },
    ),
    # from Almen 2009
    af.AnnotDef(
        name = 'ligand',
        resource = 'Almen2009',
        args = {
            'classes': 'Ligand',
        },
    ),
    af.AnnotDef(
        name = 'delta_like',
        parent = 'ligand',
        resource = 'Almen2009',
        args = {
            'classes': 'Delta',
        },
    ),
    af.AnnotDef(
        name = 'ephrin_b',
        parent = 'ligand',
        resource = 'Almen2009',
        args = {
            'classes': 'EphB',
        },
    ),
    af.AnnotDef(
        name = 'ig',
        parent = 'ligand',
        resource = 'Almen2009',
        args = {
            'classes': 'IG_Ligand',
        },
    ),
    af.AnnotDef(
        name = 'jagged',
        parent = 'ligand',
        resource = 'Almen2009',
        args = {
            'classes': 'Jagged',
        },
    ),
    af.AnnotDef(
        name = 'neuroligin',
        parent = 'ligand',
        resource = 'Almen2009',
        args = {
            'classes': 'Neuroligin',
        },
    ),
    af.AnnotDef(
        name = 'nkg2dl',
        parent = 'ligand',
        resource = 'Almen2009',
        args = {
            'classes': 'NKG2DL',
        },
    ),
    af.AnnotDef(
        name = 'semaphorin',
        parent = 'ligand',
        resource = 'Almen2009',
        args = {
            'classes': 'Semaphorins',
        },
    ),
    af.AnnotDef(
        name = 'trem_like',
        parent = 'ligand',
        resource = {'Q6UXN2'},
    ),

    # intracellular
    af.AnnotDef(
        name = 'intracellular',
        resource = af.AnnotOp(
            annots = (
                'intracellular_locate',
                'intracellular_comppi',
                'intracellular_go',
            ),
            op = set.union,
        ),
    ),
    af.AnnotDef(
        name = 'intracellular_locate',
        resource = af.AnnotOp(
            annots = (
                af.AnnotDef(
                    name = 'locate_intracellular',
                    resource = 'LOCATE',
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
                    resource = 'LOCATE',
                    args = {
                        'cls': 'cytoplasmic',
                    },
                ),
            ),
            op = set.union,
        ),
    ),
    af.AnnotDef(
        name = 'intracellular_comppi',
        resource = 'ComPPI',
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
        resource = 'GO_Intercell',
        args = {
            'mainclass': 'intracellular',
        },
    ),
    af.AnnotDef(
        name = 'intracellular_uniprot',
        resource = 'UniProt_location',
        args = {
            'location': {
                'Autophagosome',
                'Autophagosome lumen',
                'Autophagosome membrane',
                'Centriolar satellite',
                'Centriole',
                'Centromere',
                'Centrosome',
                'Chromaffin granule membrane',
                'Chromosome',
                'Cilium basal body',
                'Cis-Golgi network',
                'Cis-Golgi network membrane',
                'Cytoplasm',
                'Cytoplasmic granule',
                'Cytoplasmic granule lumen',
                'Cytoplasmic granule membrane',
                'Cytoplasmic vesicle',
                'Cytoplasmic vesicle membrane',
                'Cytoskeleton',
                'Cytosol',
                'Early endosome',
                'Early endosome membrane',
                'Endomembrane system',
                'Endoplasmic reticulum',
                'Endoplasmic reticulum lumen',
                'Endoplasmic reticulum membrane',
                'Endoplasmic reticulum-Golgi intermediate compartment', (
                    'Endoplasmic reticulum-Golgi intermediate '
                    'compartment membrane'
                ),
                'Endosome',
                'Endosome lumen',
                'Endosome membrane',
                'Golgi apparatus',
                'Golgi apparatus lumen',
                'Golgi apparatus membrane',
                'Golgi stack',
                'Golgi stack membrane',
                'Late endosome',
                'Late endosome membrane',
                'Lysosome',
                'Lysosome lumen',
                'Lysosome membrane',
                'Melanosome',
                'Melanosome membrane',
                'Microsome',
                'Microsome membrane',
                'Microtubule organizing center',
                'Mitochondrion',
                'Mitochondrion inner membrane',
                'Mitochondrion intermembrane space',
                'Mitochondrion matrix',
                'Mitochondrion membrane',
                'Mitochondrion nucleoid',
                'Mitochondrion outer membrane',
                'Multivesicular body',
                'Nuclear pore complex',
                'Nucleolus',
                'Nucleoplasm',
                'Nucleus',
                'Nucleus envelope',
                'Nucleus inner membrane',
                'Nucleus lamina',
                'Nucleus matrix',
                'Nucleus membrane',
                'Nucleus outer membrane',
                'Nucleus speckle',
                'P-body',
                'PML body',
                'Perikaryon',
                'Perinuclear region',
                'Perinuclear theca',
                'Peroxisome',
                'Peroxisome matrix',
                'Peroxisome membrane',
                'Postsynaptic Golgi apparatus',
                'Preautophagosomal structure',
                'Preautophagosomal structure membrane',
                'Recycling endosome',
                'Recycling endosome membrane',
                'Rough endoplasmic reticulum',
                'Rough endoplasmic reticulum lumen',
                'Rough endoplasmic reticulum membrane',
                'Sarcomere',
                'Sarcoplasmic reticulum',
                'Sarcoplasmic reticulum lumen',
                'Sarcoplasmic reticulum membrane',
                'Smooth endoplasmic reticulum membrane',
                'Telomere',
                'Trans-Golgi network',
            },
        },
    ),

    # extracellular
    af.AnnotDef(
        name = 'extracellular',
        resource = af.AnnotOp(
            annots = (
                'extracellular_locate',
                'extracellular_surfaceome',
                'extracellular_matrixdb',
                'extracellular_membranome',
                'extracellular_cspa',
                'extracellular_hpmr',
                'extracellular_cellphonedb',
                'extracellular_hpa',
                'extracellular_uniprot',
            ),
            op = set.union,
        ),
    ),
    af.AnnotDef(
        name = 'extracellular_locate',
        resource = af.AnnotOp(
            annots = (
                af.AnnotDef(
                    name = 'locate_extracellular',
                    resource = 'LOCATE',
                    args = {
                        'location': {
                            'extracellular',
                            'extracellular region',
                        },
                    },
                ),
                'secreted_locate',
            ),
            op = set.union,
        ),
    ),
    af.AnnotDef(
        name = 'extracellular_comppi',
        resource = 'ComPPI',
        args = {
            'location': 'extracellular',
        },
    ),
    af.AnnotDef(
        name = 'extracellular_surfaceome',
        resource = 'Surfaceome',
    ),
    af.AnnotDef(
        name = 'extracellular_matrixdb',
        resource = 'MatrixDB',
        args = {
            'mainclass': {
                'secreted',
                'ecm',
            },
        },
    ),
    af.AnnotDef(
        name = 'extracellular_cspa',
        resource = 'CSPA',
    ),
    af.AnnotDef(
        name = 'extracellular_membranome',
        resource = 'Membranome',
        args = {
            'membrane': 'Plasma membrane',
            'side': 'extracellular side',
        },
    ),
    af.AnnotDef(
        name = 'extracellular_cellphonedb',
        resource = 'CellPhoneDB',
        args = {
            'secreted': bool,
        },
    ),
    af.AnnotDef(
        name = 'extracellular_hpmr',
        resource = 'HPMR',
    ),
    af.AnnotDef(
        name = 'extracellular_membranome',
        resource = 'Membranome',
        args = {
            'membrane': 'Plasma membrane',
            'side': 'extracellular side',
        },
    ),
    af.AnnotDef(
        name = 'extracellular_hpa',
        resource = 'HPA_secretome',
        args = {
            'secreted': bool,
        },
    ),
    af.AnnotDef(
        name = 'extracellular_ramilowski',
        resource = 'Ramilowski_location',
        args = {
            'location': {
                'secreted',
            },
        },
    ),
    af.AnnotDef(
        name = 'extracellular_uniprot',
        resource = 'UniProt_location',
        args = {
            'location': {
                'Secreted',
            },
        },
    ),
    # specific subclasses from HGNC
    af.AnnotDef(
        name = 'histatin_extracellular',
        resource = 'HGNC',
        args = {
            'mainclass': 'Histatins and statherin',
        },
    ),  # secreted into saliva
    af.AnnotDef(
        name = 'proline_rich_extracellular',
        resource = 'HGNC',
        args = {
            'mainclass': 'Proline rich proteins',
        },
    ),  # secreted into saliva, enamel protective and anti-microbial
    af.AnnotDef(
        name = 'iglon5_extracellular',
        resource = {'A6NGN9'},
    ),  # function is not clear for me
    af.AnnotDef(
        name = 'immunoglobulin_heavy_extracellular',
        resource = 'HGNC',
        args = {
            'mainclass': 'Immunoglobulin heavy locus at 14q32.33',
        },
    ),  # immunoglobulin heavy chain
    af.AnnotDef(
        name = 'immunoglobulin_kappa_extracellular',
        resource = 'HGNC',
        args = {
            'mainclass': 'Immunoglobulin kappa locus at 2p11.2',
        },
    ),  # immunoglobulin V region
    af.AnnotDef(
        name = 'immunoglobulin_lambda_extracellular',
        resource = 'HGNC',
        args = {
            'mainclass': 'Immunoglobulin lambda locus at 22q11.2',
        },
    ),  # immunoglobulin V region

    # cell surface
    af.AnnotDef(
        name = 'cell_surface',
        resource = af.AnnotOp(
            annots = (
                'cell_surface_surfaceome',
                'cell_surface_go',
                'cell_surface_hpmr',
                'cell_surface_membranome',
                'cell_surface_cspa',
                'cell_surface_cellphonedb',
                'cell_surface_dgidb',
                'cell_surface_ramilowski',
                'cell_surface_uniprot',
            ),
            op = set.union,
        ),
    ),
    af.AnnotDef(
        name = 'cell_surface_hpmr',
        resource = 'HPMR',
        args = {
            'role': 'Receptor',
        },
    ),
    af.AnnotDef(
        name = 'cell_surface_surfaceome',
        resource = 'Surfaceome',
    ),
    af.AnnotDef(
        name = 'cell_surface_cspa',
        resource = 'CSPA',
    ),
    af.AnnotDef(
        name = 'cell_surface_cellphonedb',
        resource = 'CellPhoneDB',
        args = {
            'method': lambda a: a.peripheral or a.transmembrane,
        },
    ),
    af.AnnotDef(
        name = 'cell_surface_go',
        resource = 'GO_Intercell',
        args = {
            'mainclass': 'cell_surface',
        },
    ),
    af.AnnotDef(
        name = 'cell_surface_membranome',
        resource = 'Membranome',
        args = {
            'membrane': 'Plasma membrane',
            'side': 'extracellular side',
        },
    ),
    af.AnnotDef(
        name = 'cell_surface_dgidb',
        resource = 'DGIdb',
        args = {
            'category': {'CELL SURFACE', 'EXTERNAL SIDE OF PLASMA MEMBRANE'},
        },
    ),
    af.AnnotDef(
        name = 'cell_surface_ramilowski',
        resource = 'Ramilowski_location',
        args = {
            'location': 'cell surface',
        },
    ),
    af.AnnotDef(
        name = 'cell_surface_uniprot',
        resource = af.AnnotOp(
            annots = (
                af.AnnotDef(
                    name = 'cell_membrane_uniprot',
                    resource = 'UniProt_location',
                    args = {
                        'location': {
                            'Cell membrane',
                            'Acrosome membrane',
                            'Apical cell membrane',
                            'Basal cell membrane',
                            'Basolateral cell membrane',
                            'Cell surface',
                            'Dendritic spine membrane',
                            'Filopodium membrane',
                            'Flagellum membrane',
                            'Invadopodium membrane',
                            'Lamellipodium membrane',
                            'Lateral cell membrane',
                            'Microvillus membrane',
                            'Myelin membrane',
                            'Photoreceptor inner segment membrane',
                            'Presynaptic cell membrane',
                            'Ruffle membrane',
                            'Stereocilium membrane',
                        }
                    },
                ),
                # removing the transmembrane ones
                'transmembrane',
                # TODO: ensure we keep only the extracellular side
            ),
            op = set.difference,
        ),
    ),
    # specific subclasses from HGNC
    af.AnnotDef(
        name = 'glypican_cell_surface',
        resource = 'HGNC',
        args = {
            'mainclass': 'Glypicans',
        },
    ),
    af.AnnotDef(
        name = 'immunoglobulin_heavy_cell_surface',
        resource = 'HGNC',
        args = {
            'mainclass': 'Immunoglobulin heavy locus at 14q32.33',
        },
    ),  # immunoglobulin heavy chain
    af.AnnotDef(
        name = 'immunoglobulin_kappa_cell_surface',
        resource = 'HGNC',
        args = {
            'mainclass': 'Immunoglobulin kappa locus at 2p11.2',
        },
    ),  # immunoglobulin V region
    af.AnnotDef(
        name = 'immunoglobulin_lambda_cell_surface',
        resource = 'HGNC',
        args = {
            'mainclass': 'Immunoglobulin lambda locus at 22q11.2',
        },
    ),  # immunoglobulin V region
    af.AnnotDef(
        name = 'ifn_induced_cell_surface',
        resource = 'HGNC',
        args = {
            'mainclass': 'Interferon induced transmembrane proteins',
        },
    ),
    af.AnnotDef(
        name = 'ifn_induced_endocytosis_exocytosis_regulator',
        resource = 'HGNC',
        args = {
            'mainclass': 'Interferon induced transmembrane proteins',
        },
    ),

    # plasma membrane
    af.AnnotDef(
        name = 'plasma_membrane',
        resource = af.AnnotOp(
            annots = (
                'transmembrane_cellphonedb',
                'lhfpl_plasma_membrane_hgnc',
            ),
            op = set.union,
        ),
    ),
    af.AnnotDef(
        name = 'lhfpl_plasma_membrane',
        resource = 'HGNC',
        args = {
            'mainclass': 'LHFPL tetraspan proteins',
        }
    ),
    af.AnnotDef(
        name = 'plasma_membrane_regulator',
        resource = 'Almen2009',
        args = {
            'classes': 'EMP-PMP22-LIM',
        }
    ),

    # transmembrane
    af.AnnotDef(
        name = 'transmembrane',
        resource = af.AnnotOp(
            annots = (
                'transmembrane_cellphonedb',
                'transmembrane_go',
                'transmembrane_opm',
                'transmembrane_locate',
                'transmembrane_topdb',
                'transmembrane_ramilowski',
                'transmembrane_almen',
                'transmembrane_uniprot',
            ),
            op = set.union,
        ),
    ),
    af.AnnotDef(
        name = 'transmembrane',
        resource = af.AnnotOp(
            annots = (
                'transmembrane_phobius_almen',
                'transmembrane_sosui_almen',
                'transmembrane_tmhmm_almen',
            ),
            op = set.union
        ),
    ),
    af.AnnotDef(
        name = 'transmembrane_phobius',
        resource = 'Almen2009',
        args = {
            'transmembrane_phobius': True,
        },
    ),
    af.AnnotDef(
        name = 'transmembrane_sosui',
        resource = 'Almen2009',
        args = {
            'transmembrane_sosui': True,
        },
    ),
    af.AnnotDef(
        name = 'transmembrane_tmhmm',
        resource = 'Almen2009',
        args = {
            'transmembrane_tmhmm': True,
        },
    ),
    af.AnnotDef(
        name = 'transmembrane_uniprot',
        resource = af.AnnotOp(
            annots = (
                af.AnnotDef(
                    name = 'transmembrane_uniprot_location',
                    resource = 'UniProt_location',
                    args = {
                        'features': {
                            'Multi-pass membrane protein',
                            'Single-pass membrane protein',
                            'Single-pass type I membrane protein',
                            'Single-pass type II membrane protein',
                            'Single-pass type III membrane protein',
                            'Single-pass type IV membrane protein',
                        },
                    },
                ),
                af.AnnotDef(
                    name = 'transmembrane_uniprot_topology',
                    resource = 'UniProt_topology',
                    args = {
                        'topology': 'Transmembrane',
                    },
                ),
            ),
            op = set.union,
        ),
    ),
    af.AnnotDef(
        name = 'transmembrane_go',
        resource = 'GO_Intercell',
        args = {
            'mainclass': 'transmembrane',
        },
    ),
    af.AnnotDef(
        name = 'transmembrane_cellphonedb',
        resource = 'CellPhoneDB',
        args = {
            'transmembrane': True,
        },
    ),
    af.AnnotDef(
        name = 'transmembrane_opm',
        resource = 'OPM',
        args = {
            'transmembrane': True,
        },
    ),
    af.AnnotDef(
        name = 'transmembrane_topdb',
        resource = 'TopDB',
        args = {
            'topology': 'Membrane',
        },
    ),
    af.AnnotDef(
        name = 'transmembrane_locate',
        resource = 'LOCATE',
        args = {
            'cls': {
                'typeI',
                'typeII',
                'mtmp',
            },
        },
    ),
    af.AnnotDef(
        name = 'transmembrane_ramilowski',
        resource = 'Ramilowski_location',
        args = {
            'tmh': bool,
        },
    ),
    af.AnnotDef(
        name = 'lhfpl_transmembrane',
        resource = 'HGNC',
        args = {
            'mainclass': 'LHFPL tetraspan proteins',
        }
    ),

    # adhesion
    af.AnnotDef(
        name = 'adhesion',
        resource = af.AnnotOp(
            annots = (
                'adhesion_cellphonedb',
                'adhesion_go',
                'adhesion_matrisome',
                'adhesion_hgnc',
                'adhesion_integrins',
                'adhesion_zhong2015',
                'adhesion_adhesome',
                'focal_adhesion_ramilowski',
            ),
            op = set.union,
        ),
    ),
    af.AnnotDef(
        name = 'adhesion_cellphonedb',
        resource = 'CellPhoneDB',
        args = {
            'integrin': bool,
        },
    ),
    af.AnnotDef(
        name = 'adhesion_integrins',
        resource = 'Integrins',
    ),
    af.AnnotDef(
        name = 'adhesion_zhong2015',
        resource = 'Zhong2015',
    ),
    af.AnnotDef(
        name = 'adhesion',
        resource = 'HGNC',
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
        resource = 'GO_Intercell',
        args = {
            'mainclass': {'adhesion to matrix', 'adhesion to other cells'},
        }
    ),
    af.AnnotDef(
        name = 'adhesion_matrisome',
        resource = af.AnnotOp(
            annots = (
                af.AnnotDef(
                    name = 'ecm_affiliated_matrisome',
                    resource = 'Matrisome',
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
        name = 'adhesion',
        resource = 'Adhesome',
        args = {'mainclass': 'Adhesion receptor'},
    ),  # both cell-cell and cell-matrix adhesion
    af.AnnotDef(
        name = 'focal_adhesion',
        parent = 'adhesion',
        resource = 'Ramilowski_location',
        args = {
            'location': 'focal adhesion',
        },
    ),
    af.AnnotDef(
        name = 'adhesion_gprotein_coupled_receptor',
        parent = 'adhesion',
        resource = 'HGNC',
        args = {
            'mainclass': {
                'Adhesion G protein-coupled receptors, subfamily A',
                'Adhesion G protein-coupled receptors, subfamily B',
                'Adhesion G protein-coupled receptors, subfamily C',
                'Adhesion G protein-coupled receptors, subfamily D',
                'Adhesion G protein-coupled receptors, subfamily E',
                'Adhesion G protein-coupled receptors, subfamily F',
                'Adhesion G protein-coupled receptors, subfamily G',
                'Adhesion G protein-coupled receptors, subfamily L',
                'Adhesion G protein-coupled receptors, subfamily V',
            },
        }, # cell-matrix adhesion
    ),
    af.AnnotDef(
        name = '7d_cadherin',
        parent = 'cell_adhesion',
        resource = 'HGNC',
        args = {
            'mainclass': '7D cadherins',
        },
    ), # cell-cell adhesion
    af.AnnotDef(
        name = 'cadherin_related',
        parent = 'cell_adhesion',
        resource = 'HGNC',
        args = {
            'mainclass': 'Cadherin related',
        },
    ), # cell-cell adhesion
    af.AnnotDef(
        name = 'major_cadherin',
        parent = 'cell_adhesion',
        resource = 'HGNC',
        args = {
            'mainclass': 'Major cadherins',
        },
    ), # cell-cell adhesion
    af.AnnotDef(
        name = 'non_clustered_protocadherin',
        parent = 'cell_adhesion',
        resource = 'HGNC',
        args = {
            'mainclass': 'Non-clustered protocadherins',
        },
    ), # cell-cell adhesion
    af.AnnotDef(
        name = 'type1_classical_cadherin',
        parent = 'cell_adhesion',
        resource = 'HGNC',
        args = {
            'mainclass': 'Type I classical cadherins',
        },
    ), # cell-cell adhesion
    af.AnnotDef(
        name = 'type2_classical_cadherin',
        parent = 'cell_adhesion',
        resource = 'HGNC',
        args = {
            'mainclass': 'Type II classical cadherins',
        },
    ), # cell-cell adhesion
    af.AnnotDef(
        name = 'nectin',
        parent = 'cell_adhesion',
        resource = 'HGNC',
        args = {
            'mainclass': 'Nectins and nectin-like molecules',
        },
        exclude = {'O95727', 'Q15223'},
    ), # cell-cell adhesion
    af.AnnotDef(
        name = 'neurexin',
        parent = 'cell_adhesion',
        resource = 'HGNC',
        args = {
            'mainclass': 'Neurexins',
        },
    ), # cell-cell adhesion for neurons
    af.AnnotDef(
        name = 'neurexin',
        parent = 'cell_adhesion',
        resource = 'Almen2009',
        args = {
            'classes': 'Neurexin',
        },
    ),  # these are also receptors
    af.AnnotDef(
        name = 'neuroligin',
        parent = 'cell_adhesion',
        resource = 'HGNC',
        args = {
            'mainclass': 'Neuroligins',
        },
    ), # cell-cell adhesion for neurons
    af.AnnotDef(
        name = 'neuroligin',
        parent = 'cell_adhesion',
        resource = 'Almen2009',
        args = {
            'classes': 'Neuroligin',
        },
    ), # cell-cell adhesion for neurons; these are also ligands
    af.AnnotDef(
        name = 'ceacam',
        parent = 'adhesion',
        resource = 'HGNC',
        args = {
            'mainclass': (
                'Carcinoembryonic antigen related '
                'cell adhesion molecule family'
            ),
        },
    ), # in plasma membrane; cell-cell adhesion and receptors
    af.AnnotDef(
        name = 'clarin',
        parent = 'cell_adhesion',
        resource = 'HGNC',
        args = {
            'mainclass': 'Clarins',
        },
    ),  # regulation of cell-cell adhesion and synapsis in ear and retina
    af.AnnotDef(
        name = 'protocadherin',
        parent = 'cell_adhesion',
        resource = 'HGNC',
        args = {
            'mainclass': 'Clustered protocadherins',
        },
    ),  # cell-cell adhesion in brain neuronal connections
    af.AnnotDef(
        name = 'ig_like',
        parent = 'cell_adhesion',
        resource = 'HGNC',
        args = {
            'mainclass': 'Ig-like cell adhesion molecule family',
        },
        exclude = {'Q9HCN6', 'Q14CZ8'},
    ),
    af.AnnotDef(
        name = 'matrix',
        parent = 'adhesion',
        resource = {'Q9HCN6', 'Q14CZ8'},
    ),
    af.AnnotDef(
        name = 'igcam_cxadr_like',
        parent = 'cell_adhesion',
        resource = 'HGNC',
        args = {
            'mainclass': 'IgCAM CXADR-related subfamily',
        },
    ),
    af.AnnotDef(
        name = 'iglon',
        parent = 'cell_adhesion',
        resource = 'HGNC',
        args = {
            'mainclass': 'IgLON cell adhesion molecules',
        },
        exclude = {'A6NGN9'},
    ),
    af.AnnotDef(
        name = 'integrin',
        parent = 'cell_adhesion',
        resource = {
            'P23229', 'Q13349', 'Q13797', 'P20701', 'P38570', 'P05107',
            'P26010',
        },
    ),
    af.AnnotDef(
        name = 'integrin',
        parent = 'matrix_adhesion',
        resource = {
            'P06756', 'Q9UKX5', 'P08648', 'P11215', 'P26006', 'Q13683',
            'P20702', 'O75578', 'P13612', 'P17301', 'P56199', 'P08514',
            'P18564', 'O95965', 'P18084', 'P05556', 'P26012', 'P16144',
            'P05106',
        },
    ),
    af.AnnotDef(
        name = 'mucin',
        parent = 'matrix_adhesion',
        resource = af.AnnotOp(
            annots = (
                af.AnnotDef(
                    name = 'mucin_hgnc',
                    source ='HGNC',
                    args = {
                        'mainclass': 'Mucins',
                    },
                ),
                'transmembrane'
            ),
            op = set.intersection
        ),
    ),  # membrane bound mucins
    af.AnnotDef(
        name = 'receptor_tyrosine_phosphatase',
        parent = 'cell_adhesion',
        resource = {'P28827', 'O14522'},
    ),  # most of the PTPRs are not adhesion molecules but only adhesion
        # receptors or other receptors
    af.AnnotDef(
        name = 'proteoglycan',
        parent = 'cell_adhesion',
        resource = {'P16070'},
    ),
    af.AnnotDef(
        name = 'proteoglycan',
        parent = 'matrix_adhesion',
        resource = {'Q6UVK1'},
    ),
    af.AnnotDef(
        name = (
            'scavenger_receptor_cysteine_rich_'
            'matrix_adhesion_regulator_omnipath'
        ),
        resource = {'Q08380'},
    ),
    af.AnnotDef(
        name = 'scavenger_receptor_cysteine_rich',
        parent = 'cell_adhesion',
        resource = {'P30203'},
    ),
    af.AnnotDef(
        name = 'selectin',
        parent = 'cell_adhesion',
        resource = 'HGNC',
        args = {
            'mainclass': 'Selectins',
        },
    ),
    af.AnnotDef(
        name = 'adhesion',
        resource = 'Almen2009',
        args = {
            'classes': 'Adhesion',
        },
    ),
    af.AnnotDef(
        name = 'ig_like',
        parent = 'adhesion',
        resource = 'Almen2009',
        args = {
            'classes': 'IG_AdhesionProteins',
        },
    ),
    af.AnnotDef(
        name = 'mpz',
        parent = 'adhesion',
        resource = 'Almen2009',
        args = {
            'classes': 'IG_MPZ',
        },
    ),
    af.AnnotDef(
        name = 'classical_cadherin',
        parent = 'adhesion',
        resource = 'Almen2009',
        args = {
            'classes': 'CadherinClassic',
        },
    ),  # cell-cell adhesion
    af.AnnotDef(
        name = 'desmosomal_cadherin',
        parent = 'adhesion',
        resource = 'Almen2009',
        args = {
            'classes': 'CadherinOther',
        },
    ),  # cell-cell adhesion
    af.AnnotDef(
        name = 'integrin',
        parent = 'matrix_adhesion',
        resource = 'Almen2009',
        args = {
            'classes': 'Integrin',
        },
    ),  # matrix adhesion; these are also receptors
    af.AnnotDef(
        name = 'beta_protocadherin',
        parent = 'adhesion',
        resource = 'Almen2009',
        args = {
            'classes': 'ProtocadherinsBeta',
        },
    ),  # cell-cell adhesion especially between neurons
    af.AnnotDef(
        name = 'protocadherin',
        parent = 'adhesion',
        resource = 'Almen2009',
        args = {
            'classes': 'ProtocadherinsOther',
        },
    ),  # cell-cell adhesion
    af.AnnotDef(
        name = 'sarcoglycan',
        parent = 'matrix_adhesion',
        resource = 'Almen2009',
        args = {
            'classes': 'Sarcoglycan',
        },
    ),  # cell-matrix adhesion for muscle cells
    af.AnnotDef(
        name = 'selectin',
        parent = 'cell_adhesion',
        resource = 'Almen2009',
        args = {
            'classes': 'Selectin',
        },
    ),  # cell-cell adhesion between immune cells
    af.AnnotDef(
        name = 'ly6_plaur',
        parent = 'adhesion',
        resource = {'O95274', 'Q8N6Q3', 'Q8TDM5', 'Q9BY14', 'Q17RY6'},
    ),

    # surface enzyme
    af.AnnotDef(
        name = 'surface_enzyme',
        resource = af.AnnotOp(
            annots = (
                'surface_enzyme_go',
                'surface_enzyme_surfaceome',
            ),
            op = set.union,
        ),
    ),
    af.AnnotDef(
        name = 'surface_enzyme_go',
        resource = af.AnnotOp(
            annots = (
                af.AnnotOp(
                    annots = (
                        'cell_surface',
                        af.AnnotDef(
                            name = 'enzyme',
                            resource = 'GO_Intercell',
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
        resource = 'Surfaceome',
        args = {
            'mainclass': 'Enzymes',
        },
    ),
    af.AnnotDef(
        name = 'enpp_surface_enzyme',
        resource = 'HGNC',
        args = {
            'mainclass': (
                'Ectonucleotide pyrophosphatase/phosphodiesterase family'
            ),
        },
    ),  # maybe not all bound to the surface but most of them
    af.AnnotDef(
        name = 'hyaluronidase_surface_enzyme',
        resource = {'Q9UHN6', 'Q12891', 'P38567', 'Q2M3T9'},
    ),
    af.AnnotDef(
        name = 'scavenger_receptor_cysteine_rich_surface_enzyme',
        resource = {
            'P98073', 'Q9Y5Q5', 'Q9BYE2', 'Q9H3S3', 'O15393', 'P05981',
            'Q9NRS4',
        },
    ),
    af.AnnotDef(
        name = 'm1_metallopeptidase_surface_peptidase',
        resource = {'Q6Q4G3', 'Q9UIQ6', 'Q07075', 'P15144', 'Q9UKU6'},
    ),  # cleave mostly peptide ligands, hormones like TRH, angiotensin, etc
    af.AnnotDef(
        name = 'm16_metallopeptidase_surface_peptidase',
        resource = {'P14735'},
    ),  # acts on peptide hormones
    af.AnnotDef(
        name = 'm10_metallopeptidase_surface_peptidase',
        resource = {'P51511', 'P51512', 'Q9ULZ9', 'Q9Y5R2', 'Q9NPA2'},
    ),
    af.AnnotDef(
        name = 'm13_metallopeptidase_surface_peptidase',
        resource = 'HGNC',
        args = {
            'mainclass': 'M13 metallopeptidases',
        },
    ),
    af.AnnotDef(
        name = 'm14_carboxypeptidase_surface_peptidase',
        resource = {'P14384', 'O75976', 'Q8IVL8', },
    ),
    af.AnnotDef(
        name = 'vanin_surface_enzyme',
        resource = 'HGNC',
        args = {
            'mainclass': 'Vanins',
        },
        exclude = {'P43251'},
    ),

    # surface ligand
    af.AnnotDef(
        name = 'surface',
        parent = 'ligand',
        resource = af.AnnotOp(
            annots = (
                'surface_ligand_go',
                'surface_ligand_cellphonedb',
            ),
            op = set.union
        ),
        exclude = {'P0DPD6'},
    ),
    af.AnnotDef(
        name = 'surface_ligand_go',
        resource = af.AnnotOp(
            annots = (
                'cell_surface',
                'ligand_go',
            ),
            op = set.intersection,
        ),
    ),
    af.AnnotDef(
        name = 'surface_ligand_cellphonedb',
        resource = 'CellPhoneDB',
        args = {
            'method': lambda a: (
                not a.receptor and (
                    a.peripheral or
                    a.transmembrane
                )
            ),
        },
    ),
    # surface ligand subclasses from HGNC
    af.AnnotDef(
        name = 'b7_family',
        parent = 'surface_ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'B7 family',
        },
    ),
    af.AnnotDef(
        name = 'butyrophilin',
        parent = 'surface_ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Butyrophilins',
        },
    ),
    af.AnnotDef(
        name = 'ephrin',
        parent = 'surface_ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Ephrins',
        },
    ),
    af.AnnotDef(
        name = 'neuregulin',
        parent = 'surface_ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Neuregulins',
        },
    ),  # ligands for various ERBB receptors
    af.AnnotDef(
        name = 'hedgehog',
        parent = 'surface_ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Hedgehog signaling molecule family',
        },
    ),  # hedgehog proteins are initially membrane bound and can be
        # solubilized later
    af.AnnotDef(
        name = 'izumo',
        parent = 'surface_ligand',
        resource = {'Q8IYV9', 'Q6UXV1', 'Q5VZ72'},
    ),  # ligands in sperm-egg fusion
    af.AnnotDef(
        name = 'nectin',
        parent = 'surface_ligand',
        resource = {'O95727', 'Q15223'},
    ),  # ligands for T-lymphocytes
    af.AnnotDef(
        name = 'semaphorin',
        parent = 'surface_ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Semaphorins',
        },
        exclude = {
            'Q14563', 'Q13214', 'Q13275', 'Q99985',
            'Q9NS98', 'O15041', 'O95025',
        },
    ),  # surface bound ligands for plexins, regulating axonal growth
    af.AnnotDef(
        name = 'semaphorin',
        parent = 'ligand',
        resource = {
            'Q14563', 'Q13214', 'Q13275', 'Q99985',
            'Q9NS98', 'O15041', 'O95025',
        }
    ),  # secreted ligands for plexins, regulating axonal growth
    af.AnnotDef(
        name = 'anosmin',
        parent = 'ligand',
        resource = {'P23352'},
    ),
    af.AnnotDef(
        name = 'anosmin',
        parent = 'surface_ligand',
        resource = {'P23352'},
    ),
    af.AnnotDef(
        name = 'mhc',
        parent = 'surface_ligand',
        resource = 'Almen2009',
        args = {
            'classes': 'MHC',
        },
    ),
    af.AnnotDef(
        name = 'semaphorin',
        parent = 'surface_ligand',
        resource = 'Almen2009',
        args = {
            'classes': 'Semaphorin',
        },
    ),  # surface bound ligands for plexins, regulating axonal growth

    # transporter
    af.AnnotDef(
        name = 'transporter',
        resource = af.AnnotOp(
            annots = (
                'transporter_surfaceome',
                'transporter_go',
                'transporter_dgidb',
            ),
            op = set.union,
        ),
    ),
    af.AnnotDef(
        name = 'transporter_surfaceome',
        resource = 'Surfaceome',
        args = {
            'mainclass': 'Transporters',
        },
    ),
    af.AnnotDef(
        name = 'transporter_go',
        resource = 'GO_Intercell',
        args = {
            'mainclass': {
                'transport',
                'ion channels',
            },
        },
    ),
    af.AnnotDef(
        name = 'transporter_dgidb',
        resource = 'DGIdb',
        args = {
            'category': {
                'ABC TRANSPORTER',
                'TRANSPORTER',
                'ION CHANNEL',
            },
        },
    ),
    af.AnnotDef(
        name = 'transporter_adhesome',
        resource = 'Adhesome',
        args = {'mainclass': 'Channel'},
    ),
    af.AnnotDef(
        name = 'transporter',
        resource = af.AnnotOp(
            annots = (
                af.AnnotDef(
                    name = 'transporter',
                    resource = 'Almen2009',
                    args = {
                        'mainclass': 'Transporters',
                    }
                ),
                'plasma_membrane',
            ),
            op = set.intersection,
        ),
    ),
    af.AnnotDef(
        name = 'channel',
        resource = 'Adhesome',
        args = {'mainclass': 'Channel'},
    ),
    af.AnnotDef(
        name = 'channel',
        resource = af.AnnotOp(
            annots = (
                af.AnnotDef(
                    name = 'channel',
                    resource = 'Almen2009',
                    args = {
                        'classes': 'Channels',
                    }
                ),
                'plasma_membrane',
            ),
            op = set.intersection,
        ),
    ),
    af.AnnotDef(
        name = 'abca',
        parent = 'transporter',
        resource = 'Almen2009',
        args = {
            'classes': 'ABCA',
        },
    ),  # lipid transporters
    af.AnnotDef(
        name = 'abcc',
        parent = 'transporter',
        resource = 'Almen2009',
        args = {
            'classes': 'ABCC',
        },
    ),
    af.AnnotDef(
        name = 'abcg',
        parent = 'transporter',
        resource = 'Almen2009',
        args = {
            'classes': 'ABCG',
        },
    ),
    af.AnnotDef(
        name = 'solute_carrier',
        parent = 'transporter',
        resource = 'Almen2009',
        args = {
            'classes': {
                'AMAC', 'APC', 'SLC1', 'SLC10',
                'SLC11', 'SLC12', 'SLC13', 'SLC14',
                'SLC15', 'SLC16', 'SLC17', 'SLC18',
                'SLC19', 'SLC2', 'SLC20', 'SLC22',
                'SLC23', 'SLC24', 'SLC26', 'SLC27',
                'SLC28', 'SLC29', 'SLC3', 'SLC30',
                'SLC31', 'SLC34', 'SLC36', 'SLC38',
                'SLC39', 'SLC4', 'SLC40', 'SLC41',
                'SLC42(Rh)', 'SLC43', 'SLC44', 'SLC45',
                'SLC46', 'SLC5', 'SLC6', 'SLC7',
                'SLC8', 'SLC9', 'SLCO',
            },
        },
        exclude = {
            'Q8TBB6', 'Q8TE54', 'Q8IY34', 'Q9P2U8', 'Q9P2U7', 'Q8NDX2',
            'Q8NHS3', 'Q05940', 'P54219', 'Q9UGQ3', 'O95528', 'Q8N4V2',
            'Q6J4K2', 'Q71RS6', 'Q8TE54', 'Q6P1M0', 'Q6PML9', 'Q8TAD4',
            'Q6NXT4', 'O14863', 'Q8NEW0', 'Q99726', 'Q9BRI3', 'Q8NBW4',
            'Q92504', 'Q9C0K1', 'Q96H72', 'Q9UMX9', 'Q9H1V8', 'Q8TBB6',
            'Q8IVB4', 'Q9Y2E8',
        },
    ),  # transporters for various compounds, e.g. amino acids, bile acids,
        # metal ions, other inorganic and organic ions, urea, oligopeptides,
        # vitamins, sugars, steroids, organic acids, fatty acids,
        # pyrimidines, purines, nucleosides, 
        # we can split this group later
    af.AnnotDef(
        name = 'polycystin_calcium',
        parent = 'ion_channel',
        resource = 'Almen2009',
        args = {
            'classes': 'PKD1',
        },
    ),  # components of calcium channels
    af.AnnotDef(
        name = 'tmem16_calcium_dependent_chloride',
        parent = 'ion_channel',
        resource = {'Q9NQ90'},
    ),
    af.AnnotDef(
        name = 'tmem16_phospholipid_scramblase',
        resource = {'Q6IWH7', 'Q4KMQ2', 'A1A5B4'},
    ),
    af.AnnotDef(
        name = 'tmem30_aminophospholipid_flippase',
        resource = 'Almen2009',
        args = {
            'classes': 'TMEM30',
        },
    ),
    af.AnnotDef(
        name = 'tmem63_osmosensitive_cation',
        parent = 'ion_channel',
        resource = 'Almen2009',
        args = {
            'classes': 'TMEM63',
        },
        exclude = {'O94886'},
    ),

    # transporters from HGNC
    af.AnnotDef(
        name = 'water',
        parent = 'transporter',
        resource = 'HGNC',
        args = {
            'mainclass': 'Aquaporins',
        },
    ),
    af.AnnotDef(
        name = 'bestrophin',
        parent = 'transporter',
        resource = 'HGNC',
        args = {
            'mainclass': 'Bestrophins',
        },
    ),
    af.AnnotDef(
        name = 'abcc',
        parent = 'transporter',
        resource = 'HGNC',
        args = {
            'mainclass': 'ATP binding cassette subfamily C',
        },
    ),
    af.AnnotDef(
        name = 'abcg',
        parent = 'transporter',
        resource = 'HGNC',
        args = {
            'mainclass': 'ATP binding cassette subfamily G',
        },
    ),
    af.AnnotDef(
        name = 'hk_atpase',
        parent = 'transporter',
        resource = 'HGNC',
        args = {
            'mainclass': 'ATPase H+/K+ transporting',
        },
    ),
    af.AnnotDef(
        name = 'sodium_potassium_atpase',
        parent = 'transporter',
        resource = 'HGNC',
        args = {
            'mainclass': {
                'ATPase Na+/K+ transporting subunits',
                'Na+/K+ transporting ATPase interacting',
            },
        },
    ),
    af.AnnotDef(
        name = 'cnnm_metal',
        parent = 'transporter',
        resource = 'HGNC',
        args = {
            'mainclass': (
                'Cyclin and CBS domain divalent metal cation '
                'transport mediators'
            ),
        },
    ),
    af.AnnotDef(
        name = 'magnesium',
        parent = 'transporter',
        resource = 'Almen2009',
        args = {
            'classes': 'NIPA',
        },
    ),
    af.AnnotDef(
        name = 'sodium_potassium_atpase',
        parent = 'transporter',
        resource = 'Almen2009',
        args = {
            'classes': 'NKAIN',
        },
    ),
    af.AnnotDef(
        name = 'sphingolipid',
        parent = 'transporter',
        resource = 'Almen2009',
        args = {
            'classes': 'Spinster',
        },
    ),
    af.AnnotDef(
        name = 'pannexin',
        parent = 'transporter',
        resource = 'HGNC',
        args = {
            'mainclass': 'Pannexins',
        },
    ),  # as half channels they release ATP, Ca and other substances to
        # the extracellular space
    # Ion channels (subclass of transporters)
    af.AnnotDef(
        name = 'acid_sensing',
        parent = 'ion_channel',
        resource = 'HGNC',
        args = {
            'mainclass': 'Acid sensing ion channel subunits',
        },
    ),
    af.AnnotDef(
        name = 'calcium_voltage_gated',
        parent = 'ion_channel',
        resource = 'HGNC',
        args = {
            'mainclass': {
                'Calcium voltage-gated channel alpha1 subunits',
                (
                    'Calcium voltage-gated channel auxiliary '
                    'alpha2delta subunits'
                ),
                ' Calcium voltage-gated channel auxiliary beta subunits',
            },
        },
    ), # these look like all being in cell membrane, but maybe we should
       # filter them?
    af.AnnotDef(
        name = 'catsper',
        parent = 'ion_channel',
        resource = 'HGNC',
        args = {
            'mainclass': 'Cation channels sperm associated',
        },
    ),
    af.AnnotDef(
        name = 'chloride_ion_channel_regulator',
        resource = 'HGNC',
        args = {
            'mainclass': 'Chloride channel accessory',
        },
    ), # in plasma membrane, regulate cholride channels
    af.AnnotDef(
        name = 'chloride',
        parent = 'ion_channel',
        resource = 'HGNC',
        args = {
            'mainclass': 'Chloride channels, ATP-gated CFTR',
        },
    ),
    af.AnnotDef(
        name = 'cyclic_nucleotide_gated',
        parent = 'ion_channel',
        resource = 'HGNC',
        args = {
            'mainclass': 'Cyclic nucleotide gated channels',
        },
    ),  # channels gated by various cyclic nucleotides, playing roles in
        # mostly in visual and olfactory signaling
    af.AnnotDef(
        name = 'hydrogen_voltage_gated',
        parent = 'ion_channel',
        resource = 'HGNC',
        args = {
            'mainclass': 'Hydrogen voltage gated channels',
        },
    ),
    af.AnnotDef(
        name = 'b_lymphocyte_calcium',
        parent = 'ion_channel',
        resource = {'P11836'},
    ),
    af.AnnotDef(
        name = 'orai_calcium',
        parent = 'ion_channel',
        resource = 'HGNC',
        args = {
            'mainclass': 'ORAI calcium release-activated calcium modulators',
        },
    ),
    af.AnnotDef(
        name = 'calcium_activated_potassium',
        parent = 'ion_channel',
        resource = 'HGNC',
        args = {
            'mainclass': {
                'Potassium calcium-activated channel subfamily '
                'M regulatory beta subunits',
                'Potassium calcium-activated channels',
            },
        },
    ),
    af.AnnotDef(
        name = 'sodium_activated_potassium',
        parent = 'ion_channel',
        resource = 'HGNC',
        args = {
            'mainclass': 'Potassium sodium-activated channel subfamily T',
        },
    ),
    af.AnnotDef(
        name = 'two_pore_domain_potassium',
        parent = 'ion_channel',
        resource = 'HGNC',
        args = {
            'mainclass': 'Potassium two pore domain channel subfamily K',
        },
    ),
    af.AnnotDef(
        name = 'voltage_gated_potassium',
        parent = 'ion_channel',
        resource = 'HGNC',
        args = {
            'mainclass': {
                'Potassium voltage-gated channel regulatory subunits',
                'Potassium voltage-gated channel subfamily J',
                'Potassium voltage-gated channels',
            },
        },
        exclude = {
            'O43448', 'Q14722', 'Q9NZI2',
            'Q9Y2W7', 'Q13303', 'Q6PIL6',
        },
    ),
    af.AnnotDef(
        name = 'epithelial_sodium',
        parent = 'ion_channel',
        resource = 'HGNC',
        args = {
            'mainclass': 'Sodium channels epithelial',
        },
    ),
    af.AnnotDef(
        name = 'sodium_leak',
        parent = 'ion_channel',
        resource = 'HGNC',
        args = {
            'mainclass': 'Sodium leak channels, non selective',
        },
    ),
    af.AnnotDef(
        name = 'voltage_gated_sodium',
        parent = 'ion_channel',
        resource = 'HGNC',
        args = {
            'mainclass': {
                'Sodium voltage-gated channel alpha subunits',
                'Sodium voltage-gated channel beta subunits',
            },
        },
    ),
    af.AnnotDef(
        name = 'transient_receptor_potential_cation',
        parent = 'ion_channel',
        resource = 'HGNC',
        args = {
            'mainclass': 'Transient receptor potential cation channels',
        },
        exclude = {'Q9GZU1'},
    ),
    af.AnnotDef(
        name = 'transmembrane_channel_like',
        parent = 'ion_channel',
        resource = 'HGNC',
        args = {
            'mainclass': 'Transmembrane channel like family',
        },
        exclude = {'Q7Z403', 'Q8IU68'},
    ),
    af.AnnotDef(
        name = 'tweety_chloride',
        parent = 'ion_channel',
        resource = 'HGNC',
        args = {
            'mainclass': 'Tweety family',
        },
    ),  # Ca activated chloride channels
    af.AnnotDef(
        name = 'volume_regulated_anion',
        parent = 'ion_channel',
        resource = 'HGNC',
        args = {
            'mainclass': 'Volume regulated anion channel subunits',
        },
    ),
    af.AnnotDef(
        name = 'zinc_activated',
        parent = 'ion_channel',
        resource = 'HGNC',
        args = {
            'mainclass': 'Zinc activated channels',
        },
    ),
    af.AnnotDef(
        name = 'water',
        parent = 'transporter',
        resource = 'Almen2009',
        args = {
            'classes': 'Aquaporins',
        },
    ),
    af.AnnotDef(
        name = 'auxiliary',
        parent = 'transporter',
        resource = 'Almen2009',
        args = {
            'classes': 'AuxillaryTransportUnit',
        },
        exclude = {'Q9UN42'},
    ),
    af.AnnotDef(
        name = 'bestrophin',
        parent = 'transporter',
        resource = 'Almen2009',
        args = {
            'classes': 'Bestrophin',
        },
    ),
    af.AnnotDef(
        name = 'voltage_gated_calcium',
        parent = 'ion_channel',
        resource = 'Almen2009',
        args = {
            'classes': {'CACNA2D', 'CACNG'},
        },
    ),
    af.AnnotDef(
        name = 'voltage_gated_chloride',
        parent = 'ion_channel',
        resource = 'Almen2009',
        args = {
            'classes': 'CLC',
        },
    ),
    af.AnnotDef(
        name = 'calcium_activated_potassium',
        parent = 'ion_channel',
        resource = 'Almen2009',
        args = {
            'classes': 'Ca_Activated_Potassium_Channels',
        },
    ),
    af.AnnotDef(
        name = 'atp_gated',
        parent = 'ion_channel',
        resource = 'Almen2009',
        args = {
            'classes': 'ATP_gated_ion_channels',
        },
    ),
    af.AnnotDef(
        name = 'calcium',
        parent = 'ion_channel',
        resource = 'Almen2009',
        args = {
            'classes': 'CalciumChannels',
        },
    ),
    af.AnnotDef(
        name = 'inward_rectifier_potassium',
        parent = 'ion_channel',
        resource = 'Almen2009',
        args = {
            'classes': 'InwardlyRectifyingKchannel',
        },
    ),
    af.AnnotDef(
        name = 'voltage_gated_potassium',
        parent = 'ion_channel',
        resource = 'Almen2009',
        args = {
            'classes': 'KCNE',
        },
    ),
    af.AnnotDef(
        name = 'calcium_activated_potassium',
        parent = 'ion_channel',
        resource = 'Almen2009',
        args = {
            'classes': 'KCNMB',
        },
    ),
    af.AnnotDef(
        name = 'ligand_gated',
        parent = 'ion_channel',
        resource = 'Almen2009',
        args = {
            'classes': 'Ligand_gated_ion_channels',
        },
    ),
    af.AnnotDef(
        name = 'potassium',
        parent = 'ion_channel',
        resource = 'Almen2009',
        args = {
            'classes': 'Potassium_channels',
        },
    ),
    af.AnnotDef(
        name = 'beta_sodium',
        parent = 'ion_channel',
        resource = 'Almen2009',
        args = {
            'classes': 'SodiumChannelBeta',
        },
    ),
    af.AnnotDef(
        name = 'sodium',
        parent = 'ion_channel',
        resource = 'Almen2009',
        args = {
            'classes': 'SodiumChannels',
        },
    ),
    af.AnnotDef(
        name = 'transient_receptor_potential_cation',
        parent = 'ion_channel',
        resource = 'Almen2009',
        args = {
            'classes': 'TRP_channels',
        },
    ),  # calcium permeable cation channels activated by receptors
    af.AnnotDef(
        name = 'tweety_chloride',
        parent = 'ion_channel',
        resource = 'Almen2009',
        args = {
            'classes': 'Tweety',
        },
    ),
    af.AnnotDef(
        name = 'sodium_potassium',
        parent = 'ion_channel',
        resource = 'Almen2009',
        args = {
            'classes': 'Two-p_K_channels',
        },
    ),
    af.AnnotDef(
        name = 'voltage_gated',
        parent = 'ion_channel',
        resource = 'Almen2009',
        args = {
            'classes': 'Voltage_gated_ion_channels',
        },
        exclude = {'Q92736', 'P21817', 'Q14573', 'Q9ULQ1'},
    ),
    af.AnnotDef(
        name = 'xk',
        parent = 'transporter',
        resource = 'Almen2009',
        args = {
            'classes': 'XK',
        },
    ),  # transporters for amino acids,
        # scramblases for phosphatidylserine,
        # and who knows what else
    af.AnnotDef(
        name = 'sperm_associated_voltage_gated',
        parent = 'ion_channel',
        resource = 'Almen2009',
        args = {
            'classes': 'catSper-Two-P',
        },
        exclude = {'Q9ULQ1', 'Q8NHX9'},
    ),
    af.AnnotDef(
        name = 'cyclic_nucleotide_gated',
        parent = 'ion_channel',
        resource = 'Almen2009',
        args = {
            'classes': 'cycliNucleotideRegulatedChannels',
        },
    ),

    # extracellular enzyme
    af.AnnotDef(
        name = 'extracellular_enzyme',
        resource = af.AnnotOp(
            annots = (
                'extracellular_enzyme_go',
                'extracellular_enzyme_matrisome',
            ),
            op = set.union
        ),
    ),
    af.AnnotDef(
        name = 'extracellular_enzyme_go',
        resource = af.AnnotOp(
            annots = (
                # extracellular by any evidence
                # and enzyme according to GO
                af.AnnotOp(
                    annots = (
                        'extracellular',
                        af.AnnotDef(
                            name = 'enzyme',
                            resource = 'GO_Intercell',
                            args = {
                                'mainclass': 'enzyme',
                            },
                        ),
                    ),
                    op = set.intersection,
                ),
                # but not transmembrane or receptor or ligand
                af.AnnotOp(
                    annots = (
                        'transmembrane',
                        'receptor',
                        'ligand',
                    ),
                    op = set.union,
                ),
            ),
            op = set.difference,
        ),
    ),
    af.AnnotDef(
        name = 'extracellular_enzyme_matrisome',
        resource = 'Matrisome',
        args = {
            'subclass': 'ECM Regulators',
        },
    ),
    # subclasses from HGNC
    af.AnnotDef(
        name = 'adamts_extracellular_peptidase',
        resource = 'HGNC',
        args = {
            'mainclass': (
                'ADAM metallopeptidases with thrombospondin type 1 motif'
            ),
        },
    ),
    af.AnnotDef(
        name = 'adamts_like_extracellular_peptidase',
        resource = 'HGNC',
        args = {
            'mainclass': 'ADAMTS like',
        },
    ),
    af.AnnotDef(
        name = 'heparanase_extracellular_peptidase',
        resource = 'HGNC',
        args = {
            'mainclass': 'Heparanases',
        },
    ),  # act on heparin and heparane-sulphate
    af.AnnotDef(
        name = 'phospholipase',
        parent = 'extracellular_enzyme',
        resource = {
            'Q53H76', 'Q8NCC3', 'Q9NZK7', 'Q9BX93', 'P04054', 'Q13093',
            'Q5R387', 'Q9NZ20', 'P14555', 'Q9BZM1', 'Q9UNK4', 'O15496',
            'Q9BZM2', 'P39877',
        },
    ),  # secreted enzymes acting on phospholipids
    af.AnnotDef(
        name = 'defensin_extracellular_enyzme',
        resource = 'HGNC',
        args = {
            'mainclass': {
                'Defensins, alpha',
                'Defensins, beta',
            },
        },
    ),  # permeabilizing microorganism membranes or
        # binding to microorganism surfaces
    af.AnnotDef(
        name = 'lysozym_extracellular_enyzme',
        resource = 'HGNC',
        args = {
            'mainclass': {
                'Lysozymes, c-type',
                'Lysozymes, g-type'
            },
        },
    ),  # bacteriolytic proteins, some involved in sperm-egg fertilization
    af.AnnotDef(
        name = 'm14_carboxypeptidase_extracellular_peptidase',
        resource = {
            'P16870', 'P15169', 'Q8IUX7', 'Q66K79', 'P48052', 'Q96SM3',
            'Q8WXQ8', 'P15086', 'Q96IY4', 'P15085', 'Q8N4T0', 'Q9HB40',
            'P22792', 'Q9UI42', 'Q8N436', 'Q9Y646',
        },
    ),
    af.AnnotDef(
        name = 'galactosidase_extracellular_enyzme',
        resource = {'Q6UWU2', 'Q8IW92'},
    ),  # secreted galactosidases
    af.AnnotDef(
        name = 'm1_metallopeptidase_extracellular_peptidase',
        resource = {'Q9H4A4'},
    ),  #
    af.AnnotDef(
        name = 'm16_metallopeptidase_extracellular_peptidase',
        resource = {'P14735'},
    ),  # acts on peptide hormones
    af.AnnotDef(
        name = 'lipase',
        parent = 'extracellular_enzyme',
        resource = af.AnnotOp(
            annots = (
                af.AnnotDef(
                    name = 'lipase',
                    resource = 'HGNC',
                    args = {
                        'mainclass': 'Lipases',
                    },
                ),
                'extracellular',
            ),
            op = set.intersection
        ),
    ),  # secreted lipases
    af.AnnotDef(
        name = 'paraoxonase',
        parent = 'extracellular_enzyme',
        resource = 'HGNC',
        args = {
            'mainclass': 'Paraoxonases',
        },
    ),  # secreted enzymes hydrolysing lactons and other metabolites
    af.AnnotDef(
        name = 'lipocalin',
        parent = 'extracellular_enzyme',
        resource = af.AnnotOp(
            annots = (
                af.AnnotDef(
                    name = 'lipocalin',
                    resource = 'HGNC',
                    args = {
                        'mainclass': 'Lipocalins',
                    },
                ),
                'extracellular',
            ),
            op = set.intersection
        ),
    ),  # secreted lipases
    af.AnnotDef(
        name = 'immune_serin_protease_extracellular_peptidase',
        resource = 'HGNC',
        args = {
            'mainclass': (
                'Granule associated serine proteases of immune defence'
            ),
        },
    ),  # secreted by granulocytes as part of the immune response
    af.AnnotDef(
        name = 'hyaluronidase',
        parent = 'extracellular_enzyme',
        resource = {'Q12794', 'Q8WUJ3', 'O43820'},
    ),
    af.AnnotDef(
        name = (
            'inter_alpha_trypsin_inhibitor_'
            'extracellular_enyzme_regulator_hgnc'
        ),
        resource = 'HGNC',
        args = {
            'mainclass': 'Inter-alpha-trypsin inhibitor heavy chains',
        },
    ),  # protease inhibitors in plasma
    af.AnnotDef(
        name = 'm10_metallopeptidase_extracellular_peptidase',
        resource = {
            'Q9H239', 'P09237', 'P09238', 'P03956', 'P08253', 'P24347',
            'P39900', 'P45452', 'Q9NRE1', 'P22894', 'Q99542', 'Q8N119',
            'O60882', 'P14780', 'P08254',
        },
    ),  # secreted matrix metallopeptidases, many act on the ECM
    af.AnnotDef(
        name = 'ribonuclease',
        parent = 'extracellular_enzyme',
        resource = 'HGNC',
        args = {
            'mainclass': 'Ribonuclease A family',
        },
        exclude = {'P10153', 'P03950'},
    ),

    # extracellular peptidase
    af.AnnotDef(
        name = 'extracellular_peptidase_go',
        resource = af.AnnotOp(
            annots = (
                'extracellular',
                af.AnnotDef(
                    name = 'peptidase',
                    resource = 'GO_Intercell',
                    args = {
                        'mainclass': 'peptidase',
                    },
                ),
            ),
            op = set.intersection,
        ),
    ),
    af.AnnotDef(
        name = 'cela_extracellular_peptidase',
        resource = 'HGNC',
        args = {
            'mainclass': 'Chymotrypsin like elastases',
        },
    ),  # involved in ECM dynamics and remodeling
    af.AnnotDef(
        name = 'kallikrein_extracellular_peptidase',
        resource = 'HGNC',
        args = {
            'mainclass': 'Kallikreins',
        },
    ),  # extracellular serine proteases, involved in ECM dynamics
    af.AnnotDef(
        name = 'pappalysin_extracellular_peptidase',
        resource = 'HGNC',
        args = {
            'mainclass': 'Pappalysins',
        },
    ),  # cleave IGFBPs
    af.AnnotDef(
        name = (
            'scavenger_receptor_cysteine_rich_extracellular_enzyme_omnipath'
        ),
        resource = {'P58215', 'Q96JB6', 'P05156', 'P56730', 'Q96JK4', 'Q9Y4K0'},
    ),
    af.AnnotDef(
        name = (
            'glutathione_peroxidase_extracellular_enzyme_omnipath'
        ),
        resource = {'P59796', 'P49908', 'P22352'},
    ),
    af.AnnotDef(
        name = (
            'serine_peptidase_inhibitor_'
            'extracellular_peptidase_regulator_hgnc'
        ),
        resource = 'HGNC',
        args = {
            'mainclass': 'Serine peptidase inhibitors, Kazal type',
        },
    ),
    af.AnnotDef(
        name = (
            'wap_four_disulfide_core_domain_containing_'
            'extracellular_peptidase_regulator_hgnc'
        ),
        resource = 'HGNC',
        args = {
            'mainclass': 'WAP four-disulfide core domain containing',
        },
        exclude = {'P23352'},
    ),
    af.AnnotDef(
        name = (
            'tissue_inhibitor_of_metallopeptidases_'
            'extracellular_peptidase_regulator_hgnc'
        ),
        resource = 'HGNC',
        args = {
            'mainclass': 'Tissue inhibitor of metallopeptidases',
        },
    ),
    af.AnnotDef(
        name = (
            'serine_protease_extracellular_peptidase_hgnc'
        ),
        resource = 'HGNC',
        args = {
            'mainclass': 'Serine proteases',
        },
        exclude = {
            'Q9UHE8', 'P36776', 'Q9NQE7', 'Q6UWY2', 'Q7Z5A4', 'P57727',
            'Q8IU80', 'Q9Y5Q5', 'Q7RTY9', 'Q9BYE2', 'Q6ZMR5', 'Q9H3S3',
            'Q86WS5', 'Q9Y6M0', 'Q9Y5Y6', 'O60235', 'Q9NRR2', 'Q7Z410',
            'Q6ZWK6', 'Q7RTY8', 'P05981', 'Q86T26', 'Q6UWB4', 'Q9NRS4',
            'O43464', 'Q9UI38',
        },
    ),
    af.AnnotDef(
        name = (
            'serine_protease_surface_peptidase_hgnc'
        ),
        resource = {
            'Q8IU80', 'Q9Y5Q5', 'Q7RTY9', 'Q9BYE2', 'Q6ZMR5', 'Q9H3S3',
            'Q86WS5', 'Q9Y6M0', 'Q9Y5Y6', 'O60235', 'Q9NRR2', 'Q7Z410',
            'Q6ZWK6', 'Q7RTY8', 'P05981', 'Q86T26', 'Q6UWB4', 'Q9NRS4',
            'A4D1T9', 'O15393', 'Q92743', 'Q9UL52', 'Q16651'
        },
    ),
    af.AnnotDef(
        name = 'serpin_extracellular_peptidase_regulator',
        resource = af.AnnotOp(
            annots = (
                af.AnnotDef(
                    name = 'serpin_peptidase_regulator',
                    resource = 'HGNC',
                    args = {
                        'mainclass': 'Serpin peptidase inhibitors',
                    },
                ),
                'secreted',
            ),
            op = set.intersection,
        ),
    ),
    af.AnnotDef(
        name = 'biotinidase',
        parent = 'extracellular_enzyme',
        resource = {'P43251'},
    ),

    # growth factor binder or regulator
    af.AnnotDef(
        name = 'growth_factor_binder',
        resource = 'GO_Intercell',
        args = {
            'mainclass': 'growth factor binding',
        },
    ),
    af.AnnotDef(
        name = 'growth_factor_regulator',
        resource = af.AnnotOp(
            annots = (
                af.AnnotDef(
                    name = 'secreted_matrisome',
                    resource = 'Matrisome',
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
    af.AnnotDef(
        name = 'igf_binding_growth_factor_regulator',
        resource = 'HGNC',
        args = {
            'mainclass': 'Insulin like growth factor binding proteins',
        },
    ),

    # secreted
    af.AnnotDef(
        name = 'secreted',
        resource = af.AnnotOp(
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
        resource = 'LOCATE',
        args = {
            'cls': 'secretome',
        },
    ),
    af.AnnotDef(
        name = 'secreted_matrisome',
        resource = 'Matrisome',
        args = {
            'mainclass': 'Matrisome-associated',
            'subclass': 'Secreted Factors',
        },
    ),
    af.AnnotDef(
        name = 'secreted_cellphonedb',
        resource = 'CellPhoneDB',
        args = {
            'secreted': bool,
        },
    ),
    # specific subclasses from HGNC
    af.AnnotDef(
        name = 'bpi_secreted',
        resource = 'HGNC',
        args = {
            'mainclass': 'BPI fold containing',
        },
    ),

    # junctions
    # gap junction
    af.AnnotDef(
        name = 'gap_junction',
        resource = af.AnnotOp(
            annots = (
                'gap_junction_go',
                'gap_junction_ramilowski',
                'gap_junction_uniprot',
                'gap_junction_hgnc',
                'pannexin_gap_junction_hgnc',
            ),
            op = set.union,
        ),
    ),
    af.AnnotDef(
        name = 'gap_junction_go',
        resource = 'GO_Intercell',
        args = {
            'mainclass': 'gap junction',
        },
    ),
    af.AnnotDef(
        name = 'gap_junction_ramilowski',
        resource = 'Ramilowski_location',
        args = {
            'location': 'gap junction',
        },
    ),
    af.AnnotDef(
        name = 'gap_junction_uniprot',
        resource = 'UniProt_location',
        args = {
            'location': 'Gap junction',
        },
    ),
    af.AnnotDef(
        name = 'gap_junction',
        resource = 'Almen2009',
        args = {
            'classes': 'GapJunction',
        },
    ),
    af.AnnotDef(
        name = 'pannexin_gap_junction',
        resource = 'HGNC',
        args = {
            'mainclass': 'Pannexins',
        },
    ),  # either half channels or gap junctions
    # tight junction
    af.AnnotDef(
        name = 'tight_junction',
        resource = af.AnnotOp(
            annots = (
                'tight_junction_go',
                'tight_junction_ramilowski',
                'tight_junction_uniprot',
                'tight_junction_hgnc',
            ),
            op = set.union,
        ),
    ),
    af.AnnotDef(
        name = 'tight_junction_go',
        resource = 'GO_Intercell',
        args = {
            'mainclass': 'tight junction',
        },
    ),
    af.AnnotDef(
        name = 'tight_junction_ramilowski',
        resource = 'Ramilowski_location',
        args = {
            'location': 'tight junction',
        },
    ),
    af.AnnotDef(
        name = 'tight_junction_uniprot',
        resource = af.AnnotOp(
            annots = (
                af.AnnotDef(
                    name = 'tight_junction_uniprot',
                    resource = 'UniProt_location',
                    args = {
                        'location': 'Tight junction',
                    },
                ),
                'transmembrane',
            ),
            op = set.intersection
        ),
    ),
    af.AnnotDef(
        name = 'claudin_tight_junction',
        resource = 'Almen2009',
        args = {
            'mainclass': 'Claudin',
        },
    ),
    # adherens junction
    af.AnnotDef(
        name = 'adherens_junction',
        resource = af.AnnotOp(
            annots = (
                'adherens_junction_ramilowski',
                'adherens_junction_uniprot',
            ),
            op = set.union,
        ),
    ),
    af.AnnotDef(
        name = 'adherens_junction_ramilowski',
        resource = 'Ramilowski_location',
        args = {
            'location': 'adherens junction',
        },
    ),
    af.AnnotDef(
        name = 'adherens_junction_uniprot',
        resource = af.AnnotOp(
            annots = (
                af.AnnotDef(
                    name = 'adherens_junction_uniprot_location',
                    resource = 'UniProt_location',
                    args = {
                        'location': 'Adherens junction',
                    },
                ),
                af.AnnotOp(
                    annots = (
                        'transmembrane',
                        'cell_surface',
                    ),
                    op = set.union,
                ),
            ),
            op = set.intersection
        ),
    ),
    # specific subclasses from HGNC
    af.AnnotDef(
        name = 'tight_junction',
        resource = 'HGNC',
        args = {
            'mainclass': 'Claudins',
        },
    ),
    af.AnnotDef(
        name = 'desmosomal_cadherin',
        resource = 'HGNC',
        args = {
            'mainclass': 'Desmosomal cadherins',
        },
    ),
    af.AnnotDef(
        name = 'gap_junction',
        resource = 'HGNC',
        args = {
            'mainclass': 'Gap junction proteins',
        },
    ),

    # plasma membrane regions
    af.AnnotDef(
        name = 'basolateral_cell_membrane_ramilowski',
        resource = 'Ramilowski_location',
        args = {
            'location': 'basolateral cell membrane',
        },
    ),
    af.AnnotDef(
        name = 'basal_cell_membrane_ramilowski',
        resource = 'Ramilowski_location',
        args = {
            'location': 'basal cell membrane',
        },
    ),
    af.AnnotDef(
        name = 'basal_cell_membrane_ramilowski',
        resource = 'Ramilowski_location',
        args = {
            'location': 'basal cell membrane',
        },
    ),
    af.AnnotDef(
        name = 'apical_cell_membrane_ramilowski',
        resource = 'Ramilowski_location',
        args = {
            'location': 'apical cell membrane',
        },
    ),
    af.AnnotDef(
        name = 'lateral_cell_membrane_ramilowski',
        resource = 'Ramilowski_location',
        args = {
            'location': 'lateral cell membrane',
        },
    ),

    # miscellanous from HGNC
    af.AnnotDef(
        name = 'cd_molecule',
        resource = 'HGNC',
        args = {
            'mainclass': 'CD molecules',
        }, # membrane proteins, secreted, receptors, enzymes,
           # adhesion proteins
    ),
    af.AnnotDef(
        name = 'c2set_domain',
        resource = 'HGNC',
        args = {
            'mainclass': 'C2-set domain containing',
        },
    ),  # these are all plasma membrane proteins, ligands,
        # some receptors and adhesion proteins
    af.AnnotDef(
        name = 'c3_pzp_a2m',
        resource = 'HGNC',
        args = {
            'mainclass': (
                'C3 and PZP like, alpha-2-macroglobulin domain containing'
            ),
        },  # secreted or peripheral on the outer side of plasma membrane
            # enzymes, protease inhibitors, receptors, co-receptors
    ),
    af.AnnotDef(
        name = 'tetraspanin_plasma_membrane_regulator',
        resource = 'HGNC',
        args = {
            'mainclass': 'Tetraspanins',
        },
    ),  # transmembrane proteins in the plasma membrane, regulate various
        # other proteins such as channels, receptors, adhesion proteins
    af.AnnotDef(
        name = 'tetraspanin_plasma_membrane_regulator',
        resource = 'Almen2009',
        args = {
            'classes': 'Tetraspanin',
        },
        exclude = {'O60635'}
    ),  # transmembrane proteins in the plasma membrane, regulate various
        # other proteins such as channels, receptors, adhesion proteins

    # to be decided
    af.AnnotDef(
        name = 'bage',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'BAGE family',
        },
    ),
    af.AnnotDef(
        name = 'cap_ctype_lectin',
        resource = 'HGNC',
        args = {
            'mainclass': 'CAP and C-type lectin domain containing',
        },
    ), # these are secreted and affect immune signaling
    af.AnnotDef(
        name = 'cmtm',
        parent = 'receptor_regulator',
        resource = 'HGNC',
        args = {
            'mainclass': 'CKLF like MARVEL transmembrane domain containing',
        },
    ), # transmembrane in plasme membrane; regulate receptor availability
    af.AnnotDef(
        name = 'cacng_ion_channel_regulator',
        resource = 'HGNC',
        args = {
            'mainclass': ' Calcium channel auxiliary gamma subunits',
        },
    ),  # transmembrane in plasma membrane; regulate calcium channels and
        # glutamate receptors
    af.AnnotDef(
        name = 'calcium_homeostasis',
        parent = 'ion_channel',
        resource = 'HGNC',
        args = {
            'mainclass': 'Calcium homeostasis modulators',
        },
    ), # taste bud ion and ATP channels
    af.AnnotDef(
        name = 'cas_scaffold_intracell',
        parent = 'matrix_adhesion',
        resource = 'HGNC',
        args = {
            'mainclass': 'Cas scaffold proteins',
        },
    ), # intracellular part of cell-matrix (focal) adhesion signaling
    af.AnnotDef(
        name = 'cavin_caveolae',
        parent = 'intracell',
        resource = 'HGNC',
        args = {
            'mainclass': 'Cavins',
        },
    ), # caveolae formation, intercellular
    af.AnnotDef(
        name = 'clathrin_coated_pit',
        parent = 'intracell',
        resource = 'HGNC',
        args = {
            'mainclass': 'Clathrin subunits',
        },
    ), # clathrin coated pit formation, intracellular
    af.AnnotDef(
        name = 'collagen_galactosyltransferase',
        parent = 'intracell',
        resource = 'HGNC',
        args = {
            'mainclass': 'Collagen beta(1-O)galactosyltransferases',
        },
    ), # collagen synthesis (in ER)
    af.AnnotDef(
        name = 'complement_system_activator',
        resource = 'HGNC',
        args = {
            'mainclass': 'Complement system activation components',
        },
    ), # secreted receptors, enzymes and signal transmission proteins
    af.AnnotDef(
        name = 'complement_receptor_and_regulator',
        resource = 'HGNC',
        args = {
            'mainclass': 'Complement system regulators and receptors',
        },
    ),  # secreted regulators or membrane bound receptors or inhibitors
        # in the complement system downstream signaling
    af.AnnotDef(
        name = 'fibrinogen_c_domain',
        resource = 'HGNC',
        args = {
            'mainclass': 'Fibrinogen C domain containing',
        },
    ),  # all are secreted, some of them are ligands, enzymes, other kind of
        # regulators for receptors or adhesion, or ECM proteins
    af.AnnotDef(
        name = 'fibronectin_type_iii',
        resource = 'HGNC',
        args = {
            'mainclass': 'Fibronectin type III domain containing',
        },
    ),  # a mixture of plasma membrane transmembrane receptors or adhesion
        # proteins, and also ECM proteins;
        # a few of them are not extracellular at all
        # probably are annotated in other, more specific categories,
        # especially the `Ig-like cell adhesion molecule family`
    af.AnnotDef(
        name = 'immunoglobulin_like',
        resource = 'HGNC',
        args = {
            'mainclass': 'Immunoglobulin like domain containing',
        },
    ),  # a mixture of plasma membrane transmembrane receptors or adhesion
        # proteins, and also ECM proteins;
        # a few of them are not extracellular at all
        # probably are annotated in other, more specific categories,
        # especially the `Ig-like cell adhesion molecule family`
    af.AnnotDef(
        name = 'gla_domain',
        resource = 'HGNC',
        args = {
            'mainclass': 'Gla domain containing',
        },
    ),  # all secreted, various regulators of blood coagulation, ECM,
        # some enzymes or ligands or regulators of other ligands
    af.AnnotDef(
        name = 'hla',
        parent = 'surface_ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Histocompatibility complex',
        },
    ),  # histocompatibility antigen complex members for presenting
        # antigens on the cell surface
    af.AnnotDef(
        name = 'vset_domain_containing',
        resource = 'HGNC',
        args = {
            'mainclass': 'V-set domain containing',
        },
    ),  # various ligands, receptors and adhesion molecules

    # intracellular protein classes in close relation to intercellular
    # communication
    af.AnnotDef(
        name = 'crumbs_complex',
        parent = 'intracell',
        resource = 'HGNC',
        args = {
            'mainclass': 'Crumbs complex',
        },
    ),  # scaffolds and regulators for plasma membrane proteins
    af.AnnotDef(
        name = 'engulfment_motility',
        parent = 'intracell',
        resource = 'HGNC',
        args = {
            'mainclass': 'Engulfment and cell motility proteins',
        },
    ),  # some intracellular proteins involved in endocytosis
    af.AnnotDef(
        name = 'fbar_actin_dynamics_endocytosis',
        parent = 'intracell',
        resource = 'HGNC',
        args = {
            'mainclass': 'F-BAR domain containing',
        },
    ),  # intracellular proteins, most of them regulate the
        # actin dynamics in endocytosis
    af.AnnotDef(
        name = 'ferm_domain',
        parent = 'intracell',
        resource = 'HGNC',
        args = {
            'mainclass': 'FERM domain containing',
        },
    ),  # intracellular proteins, most of these regulate adhesion and
        # membrane-cytoskeleton interactions; maybe not all related closely
        # to intercellular communication processes
    af.AnnotDef(
        name = 'ferlin',
        parent = 'intracell',
        resource = 'HGNC',
        args = {
            'mainclass': 'Ferlin family',
        },
    ),  # intracellular proteins involved in plasma membrane repair
        # and synaptic vesicle fusion
    af.AnnotDef(
        name = 'fermitin',
        parent = 'intracell',
        resource = 'HGNC',
        args = {
            'mainclass': 'Fermitins',
        },
    ),  # intracellular proteins, peripheral membrane proteins on the
        # cytoplasmic side of the plasma membrane;
        # involved in cell-cell adhesion
    af.AnnotDef(
        name = 'flotillin',
        parent = 'intracell',
        resource = 'HGNC',
        args = {
            'mainclass': 'Flotillins',
        },
    ),  # intracellular proteins with a role in endocytosis
    af.AnnotDef(
        name = 'arc',
        parent = 'intracell',
        resource = {'Q7LC44'},
    ),  # intercellular RNA transfer
    af.AnnotDef(
        name = 'interferon_regulator',
        parent = 'intracell',
        resource = 'HGNC',
        args = {
            'mainclass': 'Interferon regulatory factors',
        },
    ),  # intracellular proteins mostly transcriptionally 
        # regulating interferons
    af.AnnotDef(
        name = 'junctophilin',
        parent = 'intracell',
        resource = 'HGNC',
        args = {
            'mainclass': 'Junctophilins',
        },
    ),  # intracellularily connect the plasma membrane and ER to
        # ensure quick response to membrane potential change
    af.AnnotDef(
        name = 'lims1_adhesion',
        parent = 'intracell',
        resource = {'P48059'},
    ),
    af.AnnotDef(
        name = (
            'maguk_tight_junction_intracell_omnipath'
        ),
        resource = {
            'Q07157', 'Q8N3R9', 'Q9UDY2', 'Q96QZ7', 'Q5T2T1', 'O95049',
        },
    ),  # intracellular scaffolding proteins supporting tight junctions
    af.AnnotDef(
        name = 'parin_adhesion_regulator',
        parent = 'intracell',
        resource = 'HGNC',
        args = {
            'mainclass': 'Parvins',
        },
    ),  # intracellular proteins regulating adhesion and integrin signaling
    af.AnnotDef(
        name = 'plakophilin_adhesion_regulator',
        parent = 'intracell',
        resource = 'HGNC',
        args = {
            'mainclass': 'Plakophilins',
        },
    ),  # important intracellular parts of cell-cell junctions
    af.AnnotDef(
        name = 'actin_regulation_adhesome',
        parent = 'intracell',
        resource = 'Adhesome',
        args = {'mainclass': 'Actin regulation'},
    ),
    af.AnnotDef(
        name = 'adhesion_cytoskeleton_adaptor',
        parent = 'intracell',
        resource = 'Adhesome',
        args = {'mainclass': 'Adaptor'},
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
        'interleukin_receptor_hgnc',
        'interleukin_hgnc',
        'chemokine_ligand_hgnc',
        'endogenous_ligand_hgnc',
    },
}


class_labels = {
    'ecm': 'Extracellular matrix',
    'interleukin_hgnc': 'Interleukins (HGNC)',
    'chemokine_ligand_hgnc': 'Chemokine ligands (HGNC)',
    'endogenous_ligand_hgnc': 'Endogenous ligands (HGNC)',
    'interleukin_receptor_hgnc': 'Interleukin receptors (HGNC)',
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
    'baccin': 'Baccin 2019',
    'guide2pharma': 'Guide to Pharm',
    'dgidb': 'DGIdb',
    'hpa': 'HPA',
    'lrdb': 'LRdb',
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
