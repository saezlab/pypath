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


import pypath.go as go
import pypath.dataio as dataio


"""
Gene Ontology annotations to select categories relevant in intercellular
signaling.
"""
go_annot = {
    'C': {
        'cell junction',
        'extracellular region',
        'extracellular region part',
        
        'complex of collagen trimers',
        'collagen network',
        'banded collagen fibril',
        'collagen beaded filament',
        'elastic fiber',
        'fibronectin fibril',
        ''
        
        # could not find sub-term for their membrane or lumen
        'extracellular vesicle',
        
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
        'external encapsulating structure',
        
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
    
    'M': {
        
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
        'receptor activator activity',
        'receptor inhibitor activity',
        'receptor ligand activity',
        'neurotransmitter receptor regulator activity',
        'signaling receptor activity',
        'negative regulation of signaling receptor activity',
        'positive regulation of signaling receptor activity',
        'receptor activator activity',
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
        'endocytic adaptor activity',
        
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
    
}
