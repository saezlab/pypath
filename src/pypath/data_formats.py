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

# external modules:
from future.utils import iteritems

import os
import copy

# from pypath:
import pypath.input_formats as input_formats
from pypath import common

__all__ = [
    'reaction', 'interaction', 'interaction_misc', 'pathway',
    'interaction_htp', 'ptm', 'ptm_misc', 'obsolate',
    'transcription_deprecated', 'transcription_onebyone',
    'omnipath', 'transcription', 'negative',
    'gdsc_comp_target', 'cgc', 'reactome_modifications', 'reaction_misc',
    'ligand_receptor'
]

ROOT = common.ROOT

# this is all what is needed to load the resources
# included in the pypath package
'''
Old input definitions, should not be used.
'''
obsolate = {
    'signalink2': input_formats.ReadSettings(
        name = "SignaLink2",
        separator = ",",
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = (7, ['1', '2']),
        sign = (6, '1', '-1'),
        input = os.path.join(ROOT, 'data', 'slk01human.csv'),
        references = (9, ':'),
        ncbi_tax_id = 9606,
        extra_edge_attrs = {
            "netbiol_effect": 8,
            "is_direct": 6,
            "is_directed": 7
        },
        extra_node_attrs_a = {"slk_pathways": (4, ":"),
                         "gene_name": 2},
        extra_node_attrs_b = {"slk_pathways": (5, ":"),
                         "gene_name": 3}),
    'nci_pid': input_formats.ReadSettings(
        name = "NCI-PID",
        separator = "\t",
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = (2, ['1', '-1']),
        input = os.path.join(ROOT, 'data', 'nci-pid-strict.csv'),
        references = (4, ";"),
        ncbi_tax_id = 9606,
        extra_edge_attrs = {
            "pid_effect": 2,
            "pid_evidence": (5, ";"),
            "pid_pathways": (6, ";")
        },
        extra_node_attrs_a = {},
        extra_node_attrs_b = {})
}
'''
Reaction databases.
These are not included in OmniPath, because only a minor
part of their content can be used when processing along
strict conditions to have only binary interactions with
references.
'''

reaction = {
    'Reaction resources': input_formats.ReadSettings(
        name = "reaction resources",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = (3, '1'),
        sign = None,
        resource = 4,
        input = 'get_reactions',
        references = (5, ';'),
        ncbi_tax_id = 9606,
        extra_edge_attrs = {"sif_rule": (2, ";")},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {},
        must_have_references = False)
}

pathwaycommons = {
    'PathwayCommons': input_formats.ReadSettings(
        name = "PathwayCommons",
        separator = None,
        id_col_a = 0,
        id_col_b = 2,
        id_type_a = "genesymbol",
        id_type_b = "genesymbol",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = (1, ['state-chanege', 'controls-phosphorylation-of']),
        sign = None,
        resource = 3,
        input = 'get_pathwaycommons',
        references = None,
        ncbi_tax_id = 9606,
        extra_edge_attrs = {"sif_rule": 1},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {},
        must_have_references = False,
        input_args = {
            'types': [
                'state-change', 'in-same-component', 'interacts-with',
                'controls-state-change-of', 'in-complex-with',
                'controls-transport-of', 'controls-phosphorylation-of'
            ]
        })
}

pathwaycommons1 = {
    'PathwayCommons': input_formats.ReadSettings(
        name = "PathwayCommons",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "genesymbol",
        id_type_b = "genesymbol",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = (4, '1'),
        sign = None,
        input = 'get_pathwaycommons',
        references = None,
        ncbi_tax_id = 9606,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {},
        must_have_references = False,
        input_args = {'sources_separated': False})
}

reaction_misc = {
    'nci_pid': input_formats.ReadSettings(
        name = "NCI-PID",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = (4, 'directed'),
        sign = False,
        ncbi_tax_id = 9606,
        input = 'pid_interactions',
        references = (3, ';'),
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'acsn': input_formats.ReadSettings(
        name = "ACSN",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "genesymbol",
        id_type_b = "genesymbol",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = (4, 'directed'),
        sign = False,
        ncbi_tax_id = 9606,
        input = 'acsn_interactions',
        references = (3, ';'),
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'reactome': input_formats.ReadSettings(
        name = "Reactome",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = (4, 'directed'),
        sign = False,
        ncbi_tax_id = 9606,
        huge = True,
        input = 'reactome_interactions',
        references = (3, ';'),
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {})
}

reaction_pc = {
    'acsn': input_formats.ReadSettings(
        name = "ACSN",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "genesymbol",
        id_type_b = "genesymbol",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = (2, [
            'UNKNOWN_TRANSITION', 'INTERACTION_TYPE',
            'KNOWN_TRANSITION_OMITTED', 'INHIBITION',
            'UNKNOWN_POSITIVE_INFLUENCE', 'PROTEIN_INTERACTION',
            'UNKNOWN_CATALYSIS', 'POSITIVE_INFLUENCE', 'STATE_TRANSITION',
            'TRANSLATION', 'UNKNOWN_NEGATIVE_INFLUENCE', 'NEGATIVE_INFLUENCE',
            'MODULATION', 'TRANSCRIPTION', 'COMPLEX_EXPANSION', 'TRIGGER',
            'CATALYSIS', 'PHYSICAL_STIMULATION', 'UNKNOWN_INHIBITION',
            'TRANSPORT'
        ], ';'),
        sign = (2,
              ['TRIGGER', 'UNKNOWN_POSITIVE_INFLUENCE', 'POSITIVE_INFLUENCE'],
              ['UNKNOWN_NEGATIVE_INFLUENCE', 'NEGATIVE_INFLUENCE'], ';'),
        ncbi_tax_id = 9606,
        negative_filters = [(2, ['COMPLEX_EXPANSION', 'TRANSCRIPTION'], ';'),
                         (3, 'N/A')],
        positive_filters = [],
        references = False,
        input = 'acsn_ppi',
        header = False,
        extra_edge_attrs = {'acsn_effect': (2, ';'),
                        'acsn_refs': (3, ';')},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
}
'''
Pathway databases included in OmniPath.
These are manually curated, directed, and in most
of the cases signed interactions, with literature references.
'''
pathway = {
    'trip': input_formats.ReadSettings(
        name = "TRIP",
        separator = None,
        id_col_a = 1,
        id_col_b = 0,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = (4, ['stimulation', 'inhibition']),
        sign = (4, 'stimulation', 'inhibition'),
        ncbi_tax_id = 9606,
        input = 'trip_interactions',
        references = (2, ';'),
        header = False,
        extra_edge_attrs = {'trip_methods': (3, ';')},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'spike': input_formats.ReadSettings(
        name = "SPIKE",
        separator = "\t",
        id_col_a = 1,
        id_col_b = 3,
        id_type_a = "genesymbol",
        id_type_b = "genesymbol",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = (4, ['1']),
        sign = (7, '1', '2'),
        input = 'spike_interactions',
        references = (5, ";"),
        ncbi_tax_id = 9606,
        extra_edge_attrs = {'spike_effect': 7,
                        'spike_mechanism': 11},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'signalink3': input_formats.ReadSettings(
        name = "SignaLink3",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = (6, 'directed'),
        sign = (4, 'stimulation', 'inhibition'),
        input = 'signalink_interactions',
        references = (2, ';'),
        ncbi_tax_id = 9606,
        extra_edge_attrs = {
            "netbiol_effect": 4,
            "netbiol_is_direct": 5,
            "netbiol_is_directed": 6,
            "netbiol_mechanism": 7
        },
        extra_node_attrs_a = {"slk_pathways": (8, ";")},
        extra_node_attrs_b = {"slk_pathways": (9, ";")}),
    'guide2pharma': input_formats.ReadSettings(
        name = "Guide2Pharma",
        separator = None,
        id_col_a = 0,
        id_col_b = 2,
        id_type_a = "genesymbol",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = True,
        sign = (7, 1, -1),
        input = 'get_guide2pharma',
        references = 11,
        ncbi_tax_id = 9606,
        extra_edge_attrs = {},
        extra_node_attrs_a = {'g2p_ligand_location': 8},
        extra_node_attrs_b = {'g2p_target_type': 9}),
    'ca1': input_formats.ReadSettings(
        name = "CA1",
        id_col_a = 1,
        id_col_b = 6,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = (10, ['_', '+']),
        sign = (10, '+', '_'),
        header = False,
        input = 'get_ca1',
        references = (12, ";"),
        ncbi_tax_id = 9606,
        extra_edge_attrs = {"ca1_effect": 10,
                        "ca1_type": 11},
        extra_node_attrs_a = {"ca1_location": 4,
                         "ca1_function": 3},
        extra_node_attrs_b = {"ca1_location": 9,
                         "ca1_function": 8}),
    'arn': input_formats.ReadSettings(
        name = "ARN",
        separator = ",",
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = (3, ['1', '2']),
        sign = (4, '1', '-1'),
        input = os.path.join(ROOT, 'data', 'arn_curated.csv'),
        references = (7, ":"),
        ncbi_tax_id = 9606,
        extra_edge_attrs = {
            "netbiol_effect": 4,
            "is_direct": 2,
            "is_directed": 3
        },
        extra_node_attrs_a = {"atg": 5},
        extra_node_attrs_b = {"atg": 6}),
    'nrf2': input_formats.ReadSettings(
        name = "NRF2ome",
        separator = ",",
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = (3, ['1', '2']),
        sign = (4, '1', '-1'),
        input = os.path.join(ROOT, 'data', 'nrf2ome.csv'),
        references = (5, ":"),
        ncbi_tax_id = 9606,
        extra_edge_attrs = {
            "netbiol_effect": 4,
            "is_direct": 2,
            "is_directed": 3
        },
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'macrophage': input_formats.ReadSettings(
        name = "Macrophage",
        separator = ";",
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "genesymbol",
        id_type_b = "genesymbol",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = (3, ['1']),
        sign = (2, 'Activation', 'Inhibition'),
        input = 'macrophage_interactions',
        references = (5, ","),
        ncbi_tax_id = 9606,
        extra_edge_attrs = {
            "macrophage_type": (2, ","),
            "macrophage_location": (4, ",")
        },
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'death': input_formats.ReadSettings(
        name = "DeathDomain",
        separator = "\t",
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "genesymbol",
        id_type_b = "genesymbol",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = False,
        sign = False,
        input = 'deathdomain_interactions_static',
        references = (3, ";"),
        ncbi_tax_id = 9606,
        extra_edge_attrs = {"dd_methods": (2, ';')},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'pdz': input_formats.ReadSettings(
        name = "PDZBase",
        separator = None,
        id_col_a = 1,
        id_col_b = 4,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = 1,
        sign = False,
        ncbi_tax_id = {'col': 5,
                   'dict': {
                       'human': 9606
                   }},
        input = 'get_pdzbase',
        references = 6,
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'signor': input_formats.ReadSettings(
        name = "Signor",
        separator = None,
        id_col_a = 2,
        id_col_b = 6,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        # only direct interactions
        positive_filters = [(22, 'YES')],
        # exclude TF-target interactions
        negative_filters = [(9, 'transcriptional regulation')],
        entity_type_a = "protein",
        entity_type_b = "protein",
        ncbi_tax_id = {'col': 12,
                   'dict': {
                       '9606;9606': 9606,
                       '9606': 9606
                   }},
        is_directed = (8, [
            'up-regulates', 'up-regulates activity',
            'up-regulates quantity by stabilization', 'down-regulates',
            'down-regulates activity',
            'down-regulates quantity by destabilization'
        ]),
        sign = (8, [
            'up-regulates', 'up-regulates activity',
            'up-regulates quantity by stabilization'
        ], [
            'down-regulates', 'down-regulates activity',
            'down-regulates quantity by destabilization'
        ]),
        input = 'signor_interactions',
        references = (21, ";"),
        header = True,
        extra_edge_attrs = {"signor_mechanism": (9, ';')},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {})
}

# synonym
activity_flow = pathway

'''
Interaction databases included in OmniPath.
These are subsets of the named databases, having
only low throughput, manually curated, undirected
interactions with literature references.
'''
interaction = {
    'biogrid': input_formats.ReadSettings(
        name = "BioGRID",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "genesymbol",
        id_type_b = "genesymbol",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = False,
        sign = False,
        ncbi_tax_id = 9606,
        input = 'biogrid_interactions',
        references = (2, '|'),
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'ccmap': input_formats.ReadSettings(
        name = "CancerCellMap",
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = (2, 'directed'),
        sign = False,
        ncbi_tax_id = 9606,
        input = 'get_ccmap',
        references = (3, ";"),
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'mppi': input_formats.ReadSettings(
        name = "MPPI",
        separator = "|",
        id_col_a = 2,
        id_col_b = 6,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = False,
        sign = False,
        input = 'mppi_interactions',
        references = (0, ";"),
        ncbi_tax_id = 9606,
        extra_edge_attrs = {"mppi_evidences": (1, ";")},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'dip': input_formats.ReadSettings(
        name = "DIP",
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = False,
        sign = False,
        input = 'get_dip',
        references = (2, ";"),
        ncbi_tax_id = 9606,
        extra_edge_attrs = {
            "dip_methods": (4, ";"),
            "dip_type": (3, ";"),
            'dip_id': 5
        },
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'netpath': input_formats.ReadSettings(
        name = "NetPath",
        separator = None,
        id_col_a = 1,
        id_col_b = 3,
        id_type_a = "entrez",
        id_type_b = "entrez",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = False,
        sign = False,
        input = 'netpath_interactions',
        references = (4, ";"),
        ncbi_tax_id = 9606,
        extra_edge_attrs = {
            "netpath_methods": (5, ";"),
            "netpath_type": (6, ";"),
            "netpath_pathways": (7, ';')
        },
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'innatedb': input_formats.ReadSettings(
        name = "InnateDB",
        id_col_a = 0,
        id_col_b = 2,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = False,
        sign = False,
        input = 'get_innatedb',
        references = (4, ":"),
        ncbi_tax_id = 9606,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'alz': input_formats.ReadSettings(
        name = "AlzPathway",
        separator = "\t",
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = False,
        sign = False,
        input = os.path.join(ROOT, 'data', 'alzpw-ppi.csv'),
        references = (8, ";"),
        ncbi_tax_id = 9606,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'matrixdb': input_formats.ReadSettings(
        name = "MatrixDB",
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = False,
        sign = False,
        input = 'get_matrixdb',
        references = (2, "|"),
        ncbi_tax_id = 9606,
        extra_edge_attrs = {"matrixdb_methods": (3, '|')},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
}
'''
PTM databases included in OmniPath.
These supply large sets of directed interactions.
'''
ptm = {
    'psite': input_formats.ReadSettings(
        name = "PhosphoSite",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = True,
        sign = False,
        input = 'get_phosphosite_curated',
        references = (5, ";"),
        ncbi_tax_id = 9606,
        extra_edge_attrs = {"psite_evidences": (4, ";")},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'depod': input_formats.ReadSettings(
        name = "DEPOD",
        separator = ";",
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = True,
        sign = False,
        input = 'depod_interactions',
        references = (2, "|"),
        ncbi_tax_id = 9606,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'lmpid': input_formats.ReadSettings(
        name = "LMPID",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = 0,
        sign = False,
        ncbi_tax_id = 9606,
        input = 'lmpid_interactions',
        references = (2, ';'),
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'phelm': input_formats.ReadSettings(
        name = "phosphoELM",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = 1,
        sign = False,
        ncbi_tax_id = {'col': 3,
                   'dict': {
                       'Homo sapiens': 9606
                   }},
        input = 'phelm_interactions',
        references = (2, ';'),
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'elm': input_formats.ReadSettings(
        name = "ELM",
        separator = None,
        id_col_a = 2,
        id_col_b = 3,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = 0,
        sign = False,
        ncbi_tax_id = {
            'A': {
                'col': 11,
                'dict': {
                    '"9606"(Homo sapiens)': 9606
                }
            },
            'B': {
                'col': 12,
                'dict': {
                    '"9606"(Homo sapiens)': 9606
                }
            }
        },
        input = 'get_elm_interactions',
        references = (10, ','),
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'domino': input_formats.ReadSettings(
        name = "DOMINO",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = 0,
        sign = False,
        ncbi_tax_id = {
            'A': {
                'col': 6,
                'dict': {
                    '9606': 9606
                }
            },
            'B': {
                'col': 7,
                'dict': {
                    '9606': 9606
                }
            }
        },
        input = 'domino_interactions',
        references = (5, ';'),
        header = False,
        extra_edge_attrs = {'domino_methods': (4, ';')},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'dbptm': input_formats.ReadSettings(
        name = "dbPTM",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = ['genesymbol', 'uniprot'],
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = 1,
        sign = False,
        ncbi_tax_id = 9606,
        input = 'dbptm_interactions',
        references = (2, ';'),
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {},
        must_have_references = True),
    'hprd_p': input_formats.ReadSettings(
        name = "HPRD-phos",
        separator = None,
        id_col_a = 6,
        id_col_b = 3,
        id_type_a = "genesymbol",
        id_type_b = "refseqp",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = 1,
        sign = False,
        ncbi_tax_id = 9606,
        input = 'hprd_interactions',
        references = (10, ','),
        header = False,
        extra_edge_attrs = {'hprd_mechanism': 8},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {})
}

# synonym
enzyme_substrate = ptm

'''
Other PTM datasets which are not used because the lack of
references.
'''
ptm_misc = {
    'psite_noref': input_formats.ReadSettings(
        name = "PhosphoSite_noref",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = True,
        sign = False,
        ncbi_tax_id = 9606,
        input = 'get_phosphosite_noref',
        references = False,
        extra_edge_attrs = {"psite_evidences": (4, ";")},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'ppoint': input_formats.ReadSettings(
        name = "PhosphoPoint",
        separator = ";",
        id_col_a = 1,
        id_col_b = 3,
        id_type_a = "entrez",
        id_type_b = "entrez",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = True,
        header = True,
        ncbi_tax_id = 9606,
        input = os.path.join(ROOT, 'data', 'phosphopoint.csv'),
        references = False,
        sign = False,
        extra_edge_attrs = {"phosphopoint_category": 4},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'pnetworks': input_formats.ReadSettings(
        name = "PhosphoNetworks",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "genesymbol",
        id_type_b = "genesymbol",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = True,
        sign = False,
        ncbi_tax_id = 9606,
        input = 'pnetworks_interactions',
        references = False,
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'mimp': input_formats.ReadSettings(
        name = "MIMP",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "genesymbol",
        id_type_b = "genesymbol",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = 1,
        sign = False,
        ncbi_tax_id = 9606,
        input = 'mimp_interactions',
        references = False,
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'li2012': input_formats.ReadSettings(
        name = "Li2012",
        separator = False,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "genesymbol",
        id_type_b = "genesymbol",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = 1,
        sign = False,
        input = 'li2012_interactions',
        references = False,
        ncbi_tax_id = 9606,
        extra_edge_attrs = {'li2012_mechanism': 3,
                        'li2012_route': 2},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {})
}

ptm_all = copy.deepcopy(ptm_misc)
ptm_all.update(ptm)
'''
Interaction databases not included in OmniPath.
These were omitted because lack of references,
or because we could not separate the low throughput,
manually curated interactions.
'''
interaction_misc = {
    'intact': input_formats.ReadSettings(
        name = "IntAct",
        separator = ",",
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = False,
        sign = False,
        input = 'intact_interactions',
        references = (2, ";"),
        ncbi_tax_id = 9606,
        extra_edge_attrs = {"intact_methods": (3, ';')},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'biogrid': input_formats.ReadSettings(
        name = "BioGRID",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "genesymbol",
        id_type_b = "genesymbol",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = False,
        sign = False,
        ncbi_tax_id = 9606,
        input = 'biogrid_interactions',
        references = (2, '|'),
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'hsn': input_formats.ReadSettings(
        name = "Wang",
        separator = ",",
        id_col_a = 0,
        id_col_b = 2,
        id_type_a = "entrez",
        id_type_b = "entrez",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = (4, ['Pos', 'Neg']),
        sign = (4, 'Pos', 'Neg'),
        input = 'get_hsn',
        references = False,
        ncbi_tax_id = 9606,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'acsn': input_formats.ReadSettings(
        name = "ACSN",
        separator = None,
        id_col_a = 0,
        id_col_b = 2,
        id_type_a = "genesymbol",
        id_type_b = "genesymbol",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = (1, [
            'CATALYSIS', 'UNKNOWN_CATALYSIS', 'INHIBITION',
            'PHYSICAL_STIMULATION', 'TRIGGER', 'activates',
            'UNKNOWN_POSITIVE_INFLUENCE', 'inhibits', 'MODULATION'
        ]),
        sign = (1, [
            'PHYSICAL_STIMULATION', 'TRIGGER', 'activates',
            'UNKNOWN_POSITIVE_INFLUENCE'
        ], ['INHIBITION', 'inhibits']),
        ncbi_tax_id = 9606,
        input = 'get_acsn',
        references = False,
        header = False,
        extra_edge_attrs = {'acsn_effect': 1},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'hi2': input_formats.ReadSettings(
        name = "HI-II",
        separator = None,
        id_col_a = 2,
        id_col_b = 3,
        id_type_a = "genesymbol",
        id_type_b = "genesymbol",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = False,
        sign = False,
        ncbi_tax_id = 9606,
        input = 'rolland_hi_ii_14',
        references = False,
        header = False,
        extra_edge_attrs = {'hi2_numof_screens': 4},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'hi3': input_formats.ReadSettings(
        name = "HI-III",
        separator = None,
        id_col_a = 1,
        id_col_b = 3,
        id_type_a = "genesymbol",
        id_type_b = "genesymbol",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = False,
        sign = False,
        ncbi_tax_id = 9606,
        # note: obtain the file yourself, and replace
        # this location
        input = '/home/denes/Documents/pw/data/hi3-2.3.tsv',
        references = False,
        header = True,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'lit13': input_formats.ReadSettings(
        name = "Lit-BM-13",
        separator = None,
        id_col_a = 1,
        id_col_b = 3,
        id_type_a = "genesymbol",
        id_type_b = "genesymbol",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = False,
        sign = False,
        ncbi_tax_id = 9606,
        input = 'get_lit_bm_13',
        references = False,
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'cpdb': input_formats.ReadSettings(
        name = "CPDB",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "uniprot-entry",
        id_type_b = "uniprot-entry",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = False,
        sign = False,
        ncbi_tax_id = 9606,
        input = 'get_cpdb',
        references = (3, ','),
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {})
}

interaction_htp = {
    'intact': input_formats.ReadSettings(
        name = "IntAct",
        separator = ",",
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = False,
        sign = False,
        input = 'intact_interactions',
        references = (2, ";"),
        ncbi_tax_id = 9606,
        extra_edge_attrs = {"intact_methods": (3, ';')},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {},
        input_args = {'miscore': 0.0}),
    'biogrid': input_formats.ReadSettings(
        name = "BioGRID",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "genesymbol",
        id_type_b = "genesymbol",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = False,
        sign = False,
        ncbi_tax_id = 9606,
        input = 'biogrid_interactions',
        references = (2, '|'),
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {},
        input_args = {'htp_limit': None,
                   'ltp': False}),
    'dip': input_formats.ReadSettings(
        name = "DIP",
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = False,
        sign = False,
        input = 'get_dip',
        references = (2, ";"),
        ncbi_tax_id = 9606,
        extra_edge_attrs = {
            "dip_methods": (4, ";"),
            "dip_type": (3, ";"),
            'dip_id': 5
        },
        extra_node_attrs_a = {},
        extra_node_attrs_b = {},
        input_args = {'core_only': False,
                   'small_scale_only': False}),
    'ccmap': input_formats.ReadSettings(
        name = "CancerCellMap",
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = (2, 'directed'),
        sign = False,
        ncbi_tax_id = 9606,
        input = 'get_ccmap',
        references = (3, ";"),
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'innatedb': input_formats.ReadSettings(
        name = "InnateDB",
        id_col_a = 0,
        id_col_b = 2,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = False,
        sign = False,
        input = 'get_innatedb',
        references = (4, ":"),
        ncbi_tax_id = 9606,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'matrixdb': input_formats.ReadSettings(
        name = "MatrixDB",
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = False,
        sign = False,
        input = 'get_matrixdb',
        references = (2, "|"),
        ncbi_tax_id = 9606,
        extra_edge_attrs = {"matrixdb_methods": (3, '|')},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'hprd': input_formats.ReadSettings(
        name = "HPRD",
        separator = None,
        id_col_a = 0,
        id_col_b = 3,
        id_type_a = "genesymbol",
        id_type_b = "genesymbol",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = 0,
        sign = False,
        ncbi_tax_id = 9606,
        input = 'hprd_htp',
        references = (7, ','),
        header = False,
        extra_edge_attrs = {'hprd_methods': (6, ';')},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'hi3': input_formats.ReadSettings(
        name = "Vidal HI-III",
        separator = None,
        id_col_a = 1,
        id_col_b = 3,
        id_type_a = "genesymbol",
        id_type_b = "genesymbol",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = False,
        sign = False,
        ncbi_tax_id = 9606,
        input = 'vidal_hi_iii',
        references = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {},
        must_have_references = False,
        input_args = {'fname': '/home/denes/Documents/pw/data/hi3-2.3.tsv'}),
    'mppi': input_formats.ReadSettings(
        name = "MPPI",
        separator = "|",
        id_col_a = 2,
        id_col_b = 6,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = False,
        sign = False,
        input = 'mppi_interactions',
        references = (0, ";"),
        ncbi_tax_id = 9606,
        extra_edge_attrs = {"mppi_evidences": (1, ";")},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {})
}
'''
Transcriptional regulatory interactions.
'''
transcription_onebyone = {
    'abs': input_formats.ReadSettings(
        name = "ABS",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "genesymbol",
        id_type_b = "embl_id",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = True,
        sign = False,
        ncbi_tax_id = 9606,
        input = 'get_abs',
        interaction_type = 'TF',
        references = False,
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'encode_dist': input_formats.ReadSettings(
        name = "ENCODE_distal",
        separator = None,
        id_col_a = 0,
        id_col_b = 2,
        id_type_a = "genesymbol",
        id_type_b = "genesymbol",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = True,
        sign = False,
        ncbi_tax_id = 9606,
        input = 'http://encodenets.gersteinlab.org/enets3.Distal.txt',
        interaction_type = 'TF',
        references = False,
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'encode_prox': input_formats.ReadSettings(
        name = "ENCODE_proximal",
        separator = None,
        id_col_a = 0,
        id_col_b = 2,
        id_type_a = "genesymbol",
        id_type_b = "genesymbol",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = True,
        sign = False,
        ncbi_tax_id = 9606,
        input = 'http://encodenets.gersteinlab.org/enets2.Proximal_filtered.txt',
        interaction_type = 'TF',
        references = False,
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'pazar': input_formats.ReadSettings(
        name = "PAZAR",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "enst",
        id_type_b = "ensembl",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = True,
        sign = False,
        ncbi_tax_id = 9606,
        input = 'get_pazar',
        interaction_type = 'TF',
        references = 2,
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'htri': input_formats.ReadSettings(
        name = "HTRI",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "entrez",
        id_type_b = "entrez",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = True,
        sign = False,
        ncbi_tax_id = 9606,
        input = 'get_htri',
        interaction_type = 'TF',
        references = 2,
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'oreganno': input_formats.ReadSettings(
        name = "ORegAnno",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "genesymbol",
        id_type_b = "genesymbol",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = True,
        sign = False,
        ncbi_tax_id = 9606,
        input = 'get_oreganno',
        interaction_type = 'TF',
        references = 2,
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'signor': input_formats.ReadSettings(
        name = "Signor",
        separator = None,
        id_col_a = 2,
        id_col_b = 6,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        # only direct TF-target interactions
        positive_filters = [(22, 'YES'), (9, 'transcriptional regulation')],
        entity_type_a = "protein",
        entity_type_b = "protein",
        ncbi_tax_id = {'col': 12,
                   'dict': {
                       '9606;9606': 9606,
                       '9606': 9606
                   }},
        is_directed = True,
        sign = (8, [
            'up-regulates', 'up-regulates activity',
            'up-regulates quantity by stabilization'
        ], [
            'down-regulates', 'down-regulates activity',
            'down-regulates quantity by destabilization'
        ]),
        input = 'signor_interactions',
        references = (21, ";"),
        header = True,
        extra_edge_attrs = {"signor_mechanism": (9, ';')},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {})
}

"""
New default transcription dataset is only TFregulons
as it is already an integrated resource and
has sufficient coverage.

Example
-------
import pypath

# load only `A` confidence level:
pypath.data_formats.transcription['tfregulons'].input_args['levels'] = {'A'}
pa = pypath.PyPath()
pa.init_network(pypath.data_formats.transcription)

pypath.data_formats.transcription['tfregulons'].input_args['levels'] = {
    'A', 'B', 'C', 'D'
}
pa = pypath.PyPath()
pa.init_network(pypath.data_formats.transcription)

"""
transcription = {
    'tfregulons': input_formats.ReadSettings(
        name = "TFRegulons",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "genesymbol",
        id_type_b = "genesymbol",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = True,
        sign = (2, '1', '-1'),
        ncbi_tax_id = 9606,
        input = 'get_tfregulons',
        interaction_type = 'TF',
        resource = (12, ','),
        references = (13, ','),
        header = False,
        extra_edge_attrs = {
            'tfregulons_curated': 4,
            'tfregulons_chipseq': 5,
            'tfregulons_tfbs':    6,
            'tfregulons_coexp':   7,
            'tfregulons_level':   3,
            'tfregulons_kegg_pathways': (14, '|'),
        },
        extra_node_attrs_a = {},
        extra_node_attrs_b = {})
}

'''
Old transctiptional regulation input formats.
Should not be used.
'''
transcription_deprecated = {
    'oreganno_old': input_formats.ReadSettings(
        name = "ORegAnno",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "genesymbol",
        id_type_b = "genesymbol",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = True,
        sign = False,
        ncbi_tax_id = 9606,
        input = 'get_oreganno_old',
        interaction_type = 'TF',
        references = 2,
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {})
}

'''
miRNA-target resources
'''
mirna_target = {
    'mir2dis': input_formats.ReadSettings(
        name = "miR2Disease",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "mir-mat-name",
        id_type_b = "genesymbol",
        entity_type_a = "mirna",
        entity_type_b = "protein",
        is_directed = True,
        sign = False,
        ncbi_tax_id = 9606,
        input = 'mir2disease_interactions',
        interaction_type = 'MTI',
        references = None,
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'mirdeath': input_formats.ReadSettings(
        name = "miRDeathDB",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "mirbase",
        id_type_b = "entrez",
        entity_type_a = "mirna",
        entity_type_b = "protein",
        is_directed = True,
        sign = False,
        ncbi_tax_id = {'col': 2,
                   'include': set([9606])},
        input = 'mirdeathdb_interactions',
        interaction_type = 'MTI',
        references = 3,
        header = False,
        extra_edge_attrs = {'mirdeathdb_function': 4},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'mirecords': input_formats.ReadSettings(
        name = "miRecords",
        separator = None,
        id_col_a = 0,
        id_col_b = 2,
        id_type_a = "mir-mat-name",
        id_type_b = "genesymbol",
        entity_type_a = "mirna",
        entity_type_b = "protein",
        is_directed = True,
        sign = False,
        ncbi_tax_id = {'A': {
                'col': 3,
                'dict': common.swap_dict(common.phosphoelm_taxids),
                'include': set([9606])
            },
            'B': {
                'col': 4,
                'dict': common.swap_dict(common.phosphoelm_taxids),
                'include': set([9606])
            }},
        input = 'mirecords_interactions',
        interaction_type = 'MTI',
        references = 5,
        header = False,
        extra_edge_attrs = {'mirdeathdb_function': 4},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'mirtarbase': input_formats.ReadSettings(
        name = "miRTarBase",
        separator = None,
        id_col_a = 1,
        id_col_b = 3,
        id_type_a = "mir-mat-name",
        id_type_b = "genesymbol",
        entity_type_a = "mirna",
        entity_type_b = "protein",
        is_directed = True,
        sign = False,
        ncbi_tax_id = {'A': {
                'col': 2,
                'dict': common.swap_dict(common.phosphoelm_taxids),
                'include': set([9606])
            },
            'B': {
                'col': 5,
                'dict': common.swap_dict(common.phosphoelm_taxids),
                'include': set([9606])
            }},
        positive_filters = [(7, 'Functional MTI')],
        input = 'mirtarbase_interactions',
        interaction_type = 'MTI',
        references = 8,
        header = False,
        extra_edge_attrs = {'mirtarbase_evidence': (6, '//')},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {})
}

tf_mirna = {
    'transmir': input_formats.ReadSettings(
        name = "TransmiR",
        separator = None,
        id_col_a = 0,
        id_col_b = 2,
        id_type_a = "genesymbol",
        id_type_b = "mir-mat-name",
        entity_type_a = "protein",
        entity_type_b = "mirna",
        is_directed = True,
        sign = (4, 'activation', 'repression'),
        ncbi_tax_id= {
                'col': 6,
                'dict': common.swap_dict(common.taxids),
                'include': set([9606])
            },
        input = 'transmir_interactions',
        interaction_type = 'TFMIRNA',
        references = 5,
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'encode': input_formats.ReadSettings(
        name = "ENCODE_tf-mirna",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "genesymbol",
        id_type_b = "mir-mat-name",
        entity_type_a = "protein",
        entity_type_b = "mirna",
        is_directed = True,
        sign = False,
        ncbi_tax_id = 9606,
        input = 'encode_tf_mirna_interactions',
        interaction_type = 'TFMIRNA',
        references = None,
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {})
}

lncrna_protein = {
    'lncdisease': input_formats.ReadSettings(
        name = "LncRNADisease",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "lncrna-genesymbol",
        id_type_b = "genesymbol",
        entity_type_a = "lncrna",
        entity_type_b = "protein",
        is_directed = True,
        sign = False,
        ncbi_tax_id= {
                'col': 5,
                'dict': common.swap_dict(common.taxids),
                'include': set([9606])
            },
        positive_filters = [(2, 'RNA'), (3, 'Protein')],
        input = 'lncdisease_interactions',
        interaction_type = 'LNCRP',
        references = 6,
        header = False,
        extra_edge_attrs = {'lncrnadisease_mechanism': 4},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'lncrnadb': input_formats.ReadSettings(
        name = "lncrnadb",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "lncrna-genesymbol",
        id_type_b = "genesymbol",
        entity_type_a = "lncrna",
        entity_type_b = "protein",
        is_directed = True,
        sign = False,
        ncbi_tax_id= {
                'col': 3,
                'dict': common.swap_dict(common.phosphoelm_taxids),
                'include': set([9606])
            },
        positive_filters = [(2, 'protein')],
        input = 'lncrnadb_interactions',
        interaction_type = 'LNCRP',
        references = 4,
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {})
}

ligand_receptor = {
    'ramilowski2015': input_formats.ReadSettings(
        name = "Ramilowski2015",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "genesymbol",
        id_type_b = "genesymbol",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = True,
        sign = None,
        input = 'ramilowski_interactions',
        references = (2, ','),
        ncbi_tax_id = 9606,
        extra_edge_attrs = {"ramilowski_sources": (3, ';')},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {},
        mark_source = 'ramilowski_ligand',
        mark_target = 'ramilowski_receptor',
        must_have_references = False,
        input_args = {
            'putative': False
        },
    ),
    'kirouac2010': input_formats.ReadSettings(
        name = "Kirouac2010",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "genesymbol",
        id_type_b = "genesymbol",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = True,
        sign = False,
        ncbi_tax_id = 9606,
        input = 'kirouac2010_interactions',
        references = False,
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {},
        mark_source = 'kirouac_ligand',
        mark_target = 'kirouac_receptor',
    ),
    'hpmr': input_formats.ReadSettings(
        name = 'HPMR',
        separator = None,
        id_col_a = 2,
        id_col_b = 0,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = (1, 'Ligand'),
        sign = False,
        ncbi_tax_id = 9606,
        input = 'hpmr_interactions',
        references = (3, ';'),
        must_have_references = False,
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {},
    ),
    'cellphonedb': input_formats.ReadSettings(
        name = "CellPhoneDB",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = (4, 'ligand-receptor'),
        sign = False,
        ncbi_tax_id = 9606,
        input = 'cellphonedb_interactions',
        references = (3, ';'),
        resource = (2, ';'),
        must_have_references = False,
        header = False,
        extra_edge_attrs = {'cellphonedb_type': 4},
        extra_node_attrs_a = {'cellphonedb_type': 5},
        extra_node_attrs_b = {'cellphonedb_type': 6},
        positive_filters = [],
    ),
}

small_molecule_protein = {
        'signor': input_formats.ReadSettings(
        name = "Signor",
        separator = None,
        id_col_a = 2,
        id_col_b = 6,
        id_type_a = "pubchem-cid",
        id_type_b = "uniprot",
        # only direct interactions
        positive_filters = [(22, 'YES'), (1, 'chemical')],
        # exclude TF-target interactions
        negative_filters = [(9, 'transcriptional regulation')],
        entity_type_a = "small_molecule",
        entity_type_b = "protein",
        ncbi_tax_id = {'col': 12,
                   'dict': {
                       '9606;9606': 9606,
                       '9606': 9606
                   }},
        is_directed = (8, [
            'up-regulates', 'up-regulates activity',
            'up-regulates quantity by stabilization', 'down-regulates',
            'down-regulates activity',
            'down-regulates quantity by destabilization'
        ]),
        sign = (8, [
            'up-regulates', 'up-regulates activity',
            'up-regulates quantity by stabilization'
        ], [
            'down-regulates', 'down-regulates activity',
            'down-regulates quantity by destabilization'
        ]),
        input = 'signor_interactions',
        references = (21, ";"),
        header = True,
        extra_edge_attrs = {"signor_mechanism": (9, ';')},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
}

ligand_receptor['guide2pharma'] = pathway['guide2pharma']
pathway['hpmr'] = copy.deepcopy(ligand_receptor['hpmr'])
pathway['hpmr'].must_have_references = True
pathway['hpmr'].positive_filters = []
pathway['cellphonedb'] = copy.deepcopy(ligand_receptor['cellphonedb'])
pathway['cellphonedb'].must_have_references = True
pathway['ramilowski2015'] = copy.deepcopy(ligand_receptor['ramilowski2015'])
pathway['ramilowski2015'].must_have_references = True
'''
The default set of resources in OmniPath.
'''
omnipath = {}
omnipath.update(pathway)
omnipath.update(ptm)
omnipath.update(interaction)

del omnipath['netpath']
#del omnipath['innatedb']
del omnipath['alz']
#del omnipath['biogrid']
omnipath['intact'] = interaction_htp['intact']
omnipath['biogrid'] = interaction_htp['biogrid']
omnipath['hprd'] = interaction_htp['hprd']
'''
Manually curated negative interactions, i.e. pairs of
proteins proved in experiments to not interact with
each other.
'''
negative = {
    'negatome': input_formats.ReadSettings(
        name = "Negatome",
        separator = "\t",
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = 0,
        input = 'negatome_pairs',
        ncbi_tax_id = 9606,
        extra_edge_attrs = {"references": (2, ';'),
                        "negatome_methods": (3, ';')},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {})
}

biocarta = input_formats.ReadSettings(
    name = "BioCarta",
    separator = ";",
    id_col_a = 0,
    id_col_b = 2,
    id_type_a = "entrez",
    id_type_b = "entrez",
    entity_type_a = "protein",
    entity_type_b = "protein",
    is_directed = 1,
    input = os.path.join(ROOT, 'data', 'biocarta-pid.csv'),
    extra_edge_attrs = {},
    extra_node_attrs_a = {"biocarta_pathways": (4, ",")},
    extra_node_attrs_b = {"biocarta_pathways": (4, ",")})

nci_pid = input_formats.ReadSettings(
    name = "NCI-PID",
    separator = ";",
    id_col_a = 0,
    id_col_b = 1,
    id_type_a = "uniprot",
    id_type_b = "uniprot",
    entity_type_a = "protein",
    entity_type_b = "protein",
    is_directed = 1,
    input = os.path.join(ROOT, 'data', 'nci-pid.csv'),
    extra_edge_attrs = {},
    extra_node_attrs_a = {"nci_pid_pathways": (2, ",")},
    extra_node_attrs_b = {"nci_pid_pathways": (2, ",")})

reactome = input_formats.ReadSettings(
    name = "Reactome",
    separator = ";",
    id_col_a = 0,
    id_col_b = 1,
    id_type_a = "uniprot",
    id_type_b = "uniprot",
    entity_type_a = "protein",
    entity_type_b = "protein",
    is_directed = 1,
    input = os.path.join(ROOT, 'data', 'reactome-pid.csv'),
    extra_edge_attrs = {},
    extra_node_attrs_a = {"reactome_pathways": (2, ",")},
    extra_node_attrs_b = {"reactome_pathways": (2, ",")})

gdsc_comp_target = input_formats.ReadSettings(
    name = "GDSC",
    separator = ";",
    id_col_a = 1,
    id_col_b = 0,
    id_type_a = "pubchem",
    id_type_b = "genesymbol",
    entity_type_a = "drug",
    entity_type_b = "protein",
    is_directed = 1,
    input = "gdsc.sif",
    extra_edge_attrs = {},
    extra_node_attrs_a = {"gene_name": 2},
    extra_node_attrs_b = {})

gdsc_lst = input_formats.ReadList(
    name = "GDSC",
    separator = ";",
    id_col = 0,
    id_type = "genesymbol",
    entity_type = "protein",
    input = os.path.join(ROOT, 'data', 'gdsc.sif'),
    extra_attrs = {'drugs': 2})

gdsc_lst = input_formats.ReadList(
    name = "atg",
    separator = ";",
    id_col = 0,
    id_type = "genesymbol",
    entity_type = "protein",
    input = os.path.join(ROOT, 'data', 'autophagy.list'),
    extra_attrs = {'drugs': 2})

cgc = input_formats.ReadList(
    name = "cgc",
    id_col = 2,
    id_type = "entrez",
    entity_type = "protein",
    input = 'get_cgc',
    extra_attrs = {})

intogen_cancer = input_formats.ReadList(
    name = "IntOGen",
    separator = "\t",
    id_col = 1,
    id_type = "genesymbol",
    entity_type = "protein",
    input = None,
    extra_attrs = {})

reactome_modifications = {
    'phosphorylated': ('phosphorylation', 'X'),
    'glycosylated': ('glycosylation', 'X'),
    'acetylated': ('acetylated', 'X'),
    'prenylated': ('prenylation', 'X'),
    'ubiquitinated': ('ubiquitination', 'X'),
    'myristoylated': ('myristoylation', 'X'),
    'hydroxylated': ('hydroxylation', 'X'),
    'acetylated residue': ('acetylation', 'X'),
    'palmitoylated residue': ('palmitoylation', 'X'),
    'sumoylated lysine': ('sumoylation', 'K'),
    'O-palmitoyl-L-threonine': ('palmitoylation', 'T'),
    'acetylated L-serine': ('acetylation', 'S'),
    'glycosylated residue': ('glycosylation', 'X'),
    'methylated L-arginine': ('methylation', 'R'),
    'ubiquitination': ('ubiquitination', 'X'),
    'phosphorylated residue': ('phosphorylation', 'X'),
    'O-phospho-L-threonine': ('phosphorylation', 'T'),
    'O-glycosyl-L-threonine': ('glycosylation', 'T'),
    'methylated L-lysine': ('methylation', 'K'),
    'myristoylated residue': ('myristoylation', 'X'),
    'N-myristoyl-glycine': ('myristoylation', 'G'),
    'O-palmitoyl-L-serine': ('palmitoylation', 'S'),
    'palmitoylated residue [residue = N]': ('palmitoylation', 'X'),
    'N-acetylated L-lysine': ('acetylation', 'K'),
    'O-glycosyl-L-serine': ('glycosylation', 'S'),
    'N-acetyl-L-methionine': ('acetylation', 'M'),
    'ubiquitinylated lysine': ('ubiquitination', 'K'),
    'S-farnesyl-L-cysteine': ('farnesylation', 'C'),
    'S-phospho-L-cysteine': ('phosphorylation', 'C'),
    'hydroxylated proline': ('hydroxylation', 'P'),
    'palmitoylated residue [residue = Y]': ('palmitoylation', 'Y'),
    'O4\'-phospho-L-tyrosine': ('phosphorylation', 'Y'),
    'O-phospho-L-serine': ('phosphorylation', 'S'),
    'O-phospho-L-threonine': ('phosphorylation', 'T'),
    '(2S,4R)-4-hydroxyproline': ('hydroxylation', 'P'),
    '(2S,3S)-3-hydroxyproline': ('hydroxylation', 'P'),
    'O5-galactosyl-L-hydroxylysine': ('galactosytlation', 'K'),
    '(2S,5R)-5-hydroxylysine': ('hydroxylation', 'K'),
    'O5-glucosylgalactosyl-L-hydroxylysine': ('glucosylgalactosylation', 'K'),
    'N4-glycosyl-L-asparagine': ('glycosylation', 'N'),
    'S-palmitoyl-L-cysteine': ('palmitoylation', 'C'),
    'N-myristoylglycine': ('myristoylation', 'G'),
    'half cystine': ('half cystine', 'C'),
    'S-geranylgeranyl-L-cysteine': ('geranylation', 'C'),
    'N6-acetyl-L-lysine': ('acetylation', 'K'),
    'N\'-formyl-L-kynurenine': ('formylation', 'W'),
    'Oxohistidine (from histidine)': ('oxo', 'H'),
    'dihydroxyphenylalanine (Phe)': ('dihydroxylation', 'F'),
    'glutamyl semialdehyde (Pro)': ('glutamylation', 'P'),
    'monohydroxylated asparagine': ('hydroxylation', 'N'),
    'monohydroxylated proline': ('hydroxylation', 'P'),
    'ubiquitinylated lysine': ('ubiquitination', 'K'),
    'N6,N6,N6-trimethyl-L-lysine': ('trimethylation', 'K'),
    'N6,N6-dimethyl-L-lysine': ('dimethylation', 'K'),
    'N6-myristoyl-L-lysine': ('myristoylation', 'K'),
    'sumoylated lysine': ('sumoylation', 'K'),
    'N6-methyl-L-lysine': ('methylation', 'K'),
    'omega-N-methyl-L-arginine': ('methylation', 'R'),
    'asymmetric dimethyl-L-arginine': ('dimethylation', 'R'),
    'symmetric dimethyl-L-arginine': ('dimethylation', 'R'),
    'O4\'-glucosyl-L-tyrosine': ('glycosylation', 'Y'),
    'N6-biotinyl-L-lysine': ('biotinylation', 'K'),
    'O-acetyl-L-serine': ('acetylation', 'S'),
    '1-thioglycine': ('thiolation', 'G'),
    'S-acetyl-L-cysteine': ('acetylation', 'C'),
    'N-acetyl-L-alanine': ('acetylation', 'A'),
    'S-methyl-L-cysteine': ('methylation', 'C'),
    'L-gamma-carboxyglutamic acid': ('carboxylation', 'Z'),
    '(2S,3R)-3-hydroxyaspartic acid': ('hydroxylation', 'D'),
    'O-fucosyl-L-threonine': ('fucosylation', 'T'),
    'O-fucosyl-L-serine': ('fucosylation', 'S'),
    'O-palmitoleyl-L-serine': ('palmitoylation', 'S'),
    '1-thioglycine (C-terminal)': ('thiolation', 'G'),
    'neddylated lysine': ('neddylation', 'K'),
    'N-palmitoyl-L-cysteine': ('palmitoylation', 'C'),
    'S-farnesyl-L-cysteine': ('farnesylation', 'C')
}

categories = {
    'Vidal HI-III': 'i',
    'CancerCellMap': 'p',
    'InnateDB': 'i',
    'SPIKE': 'p',
    'LMPID': 'm',
    'DIP': 'i',
    'HPRD': 'i',
    'HPRD-phos': 'm',
    'PDZBase': 'p',
    'dbPTM': 'm',
    'MatrixDB': 'i',
    'DOMINO': 'm',
    'Signor': 'p',
    'Macrophage': 'p',
    'NetPath': 'r',
    'ELM': 'm',
    'SignaLink2': 'p',
    'SignaLink3': 'p',
    'NRF2ome': 'p',
    'DEPOD': 'm',
    'phosphoELM': 'm',
    'MPPI': 'i',
    'Guide2Pharma': 'l',
    'Guide2Pharma_CP': 'l',
    'TRIP': 'p',
    'AlzPathway': 'r',
    'PhosphoSite': 'm',
    'CA1': 'p',
    'NCI-PID': 'r',
    'DeathDomain': 'p',
    'ARN': 'p',
    'BioGRID': 'i',
    'IntAct': 'i',
    'Reactome': 'r',
    'ACSN': 'r',
    'WikiPathways': 'r',
    'PANTHER': 'r',
    'ABS': 't',
    'ENCODE_distal': 't',
    'PAZAR': 't',
    'ENCODE_proximal': 't',
    'ORegAnno': 't',
    'HTRI': 't',
    'MIMP': 'm',
    'PhosphoNetworks': 'm',
    'Li2012': 'm',
    'PhosphoPoint': 'm',
    'PhosphoSite_noref': 'm',
    'Ramilowski2015': 'l',
    'Kirouac2010': 'l',
    'HPMR': 'l',
    'CellPhoneDB': 'l',
    'Guide2Pharma': 'l',
    'GO_lig_rec': 'l',
    'guidetopharmacology.org': 'l',
    'UniProt': 'l',
    'InnateDB-All': 'i',
    'MINT': 'i',
}

p = set([])
i = set([])
r = set([])
m = set([])
t = set([])
l = set([])

for db, cats in iteritems(categories):
    
    for c in cats:
        
        locals()[c].add(db)

catnames = {
    'm': 'Enzyme-substrate',
    'p': 'Activity flow',
    'i': 'Undirected PPI',
    'r': 'Process description',
    't': 'Transcription',
    'l': 'Ligand-receptor',
}

catletters = dict(map(reversed, iteritems(catnames)))

pathway_resources = p
interaction_resources = i
ptm_resources = m
reaction_resources = r
transctiption_resources = t
ligand_receptor_resources = l
