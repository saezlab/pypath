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
from past.builtins import xrange, range

import imp

import pypath.go as go
import pypath.dataio as dataio


"""
Gene Ontology annotations to select categories relevant in intercellular
signaling.
"""

intercell_categories = {
    'junction': {
        'C': {'cell_junction'},
    },
}

intercell_go_terms = {
    
    # cellular component
    'C': {
        'junction': {
            'cell junction',
        },
        'extracellular': {
            'extracellular region',
            'extracellular region part',
        },
        'extracellular_matrix': {
            'extracellular matrix',
            'complex of collagen trimers',
            'collagen network',
            'banded collagen fibril',
            'collagen beaded filament',
            'elastic fiber',
            'fibronectin fibril',
        },
        'exosome': {
            # could not find sub-term for their membrane or lumen
            'extracellular vesicle',
        },
        'cell_surface': {
            # only plasma membrane components facing outside
            'cell surface',
            'external side of plasma membrane',
        },
        ''
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


class Intercell(object):
    
    def __init__(
            self,
            pa = None,
            annot = None,
            categories = None,
        ):
        
        self.pa = pa
        self.annot = (
            annot
                if annot else
            pa.go[pa.ncbi_tax_id]
                if hasattr(pa, 'go') else
            go.GOAnnotation()
        )
        
        self.names = categories or categories_default
    
    def reload(self):
        """
        Reloads the object from the module level.
        """
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    def get_go_ids(self):
        
        self.terms = dict(
            (
                domain,
                set(self.annot.ontology.get_term(name) for name in names)
            )
            for domain, names in iteritems(self.names)
        )
