#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module.
#  Contains meta-information on all databases/resources.
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

#TODO will be moved to jsons

# external modules:
from future.utils import iteritems

import os
import copy

# from pypath:
import pypath.internals.input_formats as input_formats
import pypath.resources.urls as urls
import pypath.share.common as common
import pypath.utils.taxonomy as taxonomy


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
    'signalink2': input_formats.NetworkInput(
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
    'nci_pid': input_formats.NetworkInput(
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
    'Reaction resources': input_formats.NetworkInput(
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

pathwaycommons_all = input_formats.NetworkInput(
    name = "PathwayCommons",
    separator = None,
    id_col_a = 0,
    id_col_b = 2,
    id_type_a = "genesymbol",
    id_type_b = "genesymbol",
    entity_type_a = "protein",
    entity_type_b = "protein",
    is_directed = (
        1,
        [
            'state-chanege',
            'controls-phosphorylation-of',
            'controls-state-change-of',
            'controls-transport-of',
        ],
    ),
    sign = None,
    input = 'pathwaycommons.pathwaycommons_interactions',
    references = None,
    ncbi_tax_id = 9606,
    resource = 3,
    extra_edge_attrs = {"pc_rule": 1},
    extra_node_attrs_a = {},
    extra_node_attrs_b = {},
    must_have_references = False,
    input_args = {
        'types': {
            'state-change',
            'in-same-component',
            'interacts-with',
            'controls-state-change-of',
            'in-complex-with',
            'controls-transport-of',
            'controls-phosphorylation-of',
        },
    }
)


def _pathwaycommons_single_resource(resource):

    dmodel_interaction = {'CORUM', 'IntAct', 'DIP', 'BioGRID', 'BIND', 'INOH'}

    input_def = copy.deepcopy(pathwaycommons_all)
    input_def.input_args['resources'] = resource
    input_def.data_model = (
        'interaction'
            if resource in dmodel_interaction else
        'activity_flow'
    )

    return input_def


pathwaycommons = dict(
    (
        resource.lower().replace('-', '_'),
        _pathwaycommons_single_resource(resource),

    )
    for resource in (
        'NCI-PID',
        'KEGG',
        'CORUM',
        'BIND',
        'HPRD',
        'WikiPathways',
        'INOH',
        'BioGRID',
        'NetPath',
        'Reactome',
        'DIP',
        'IntAct',
        'PANTHER',
        'PhosphoSite',
    )
)

# synonym for old name
reactome_pc = pathwaycommons['reactome']


reaction_misc = {
    'nci_pid': input_formats.NetworkInput(
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
    'acsn': input_formats.NetworkInput(
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
    'reactome': input_formats.NetworkInput(
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
        input_args = {'ask': False},
        references = (3, ';'),
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {})
}

reaction_pc = {
    'acsn': input_formats.NetworkInput(
        name = 'ACSN',
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = 'genesymbol',
        id_type_b = 'genesymbol',
        entity_type_a = 'protein',
        entity_type_b = 'protein',
        is_directed = (
            2,
            [
                'UNKNOWN_TRANSITION',
                'INTERACTION_TYPE',
                'KNOWN_TRANSITION_OMITTED',
                'INHIBITION',
                'UNKNOWN_POSITIVE_INFLUENCE',
                'UNKNOWN_CATALYSIS',
                'POSITIVE_INFLUENCE',
                'STATE_TRANSITION',
                'TRANSLATION',
                'UNKNOWN_NEGATIVE_INFLUENCE',
                'NEGATIVE_INFLUENCE',
                'MODULATION',
                'TRANSCRIPTION',
                'COMPLEX_EXPANSION',
                'TRIGGER',
                'CATALYSIS',
                'PHYSICAL_STIMULATION',
                'UNKNOWN_INHIBITION',
                'TRANSPORT',
                'inhibits',
                'activates',
            ],
            ';',
        ),
        sign = (
            2,
            [
                'TRIGGER',
                'UNKNOWN_POSITIVE_INFLUENCE',
                'POSITIVE_INFLUENCE',
                'PHYSICAL_STIMULATION',
                'activates',
            ],
            [
                'INHIBITION',
                'UNKNOWN_NEGATIVE_INFLUENCE',
                'NEGATIVE_INFLUENCE',
                'inhibits',
            ],
            ';',
        ),
        ncbi_tax_id = 9606,
        negative_filters = [
            (
                2,
                [
                    'COMPLEX_EXPANSION',
                    'TRANSCRIPTION',
                ],
                ';'
            ),
        ],
        positive_filters = [],
        references = (3, ';'),
        input = 'acsn_interactions',
        header = False,
        extra_edge_attrs = {'acsn_effect': (2, ';')},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
}


'''
Pathway databases included in OmniPath.
These are manually curated, directed, and in most
of the cases signed interactions, with literature references.
'''
pathway = {
    'trip': input_formats.NetworkInput(
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
        input = 'trip.trip_interactions',
        references = (2, ';'),
        header = False,
        extra_edge_attrs = {'trip_methods': (3, ';')},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}
    ),
    'spike': input_formats.NetworkInput(
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
        input = 'spike.spike_interactions',
        references = (5, ";"),
        ncbi_tax_id = 9606,
        extra_edge_attrs = {
            'spike_effect': 7,
            'spike_mechanism': 11,
        },
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}
    ),
    'signalink3': input_formats.NetworkInput(
        name = 'SignaLink3',
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = 'uniprot',
        id_type_b = 'uniprot',
        entity_type_a = 'protein',
        entity_type_b = 'protein',
        is_directed = (3, True),
        sign = (4, 1, -1),
        input = 'signalink.signalink_interactions',
        references = 9,
        resource = 10,
        ncbi_tax_id = 9606,
        extra_edge_attrs = {},
        positive_filters = [(2, True)], # only direct interactions
        extra_node_attrs_a = {},
        extra_node_attrs_b = {},
    ),
    'guide2pharma': input_formats.NetworkInput(
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
        input = 'guide2pharma_interactions',
        references = 11,
        ncbi_tax_id = 9606,
        extra_edge_attrs = {},
        extra_node_attrs_a = {'g2p_ligand_location': 8},
        extra_node_attrs_b = {'g2p_target_type': 9},
    ),
    'ca1': input_formats.NetworkInput(
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
    'arn': input_formats.NetworkInput(
        name = "ARN",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = (3, ['1', '2']),
        sign = (4, '1', '-1'),
        input = 'arn_interactions',
        references = (7, ":"),
        ncbi_tax_id = 9606,
        extra_edge_attrs = {
            "netbiol_effect": 4,
            "is_direct": 2,
            "is_directed": 3
        },
        extra_node_attrs_a = {"atg": 5},
        extra_node_attrs_b = {"atg": 6}),
    'nrf2': input_formats.NetworkInput(
        name = "NRF2ome",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = (3, ['1', '2']),
        sign = (4, '1', '-1'),
        input = 'nrf2ome_interactions',
        references = (5, ":"),
        ncbi_tax_id = 9606,
        extra_edge_attrs = {
            "netbiol_effect": 4,
            "is_direct": 2,
            "is_directed": 3
        },
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'macrophage': input_formats.NetworkInput(
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
    'death': input_formats.NetworkInput(
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
        input = 'deathdomain.deathdomain_interactions_rescued',
        references = (3, ";"),
        ncbi_tax_id = 9606,
        extra_edge_attrs = {"dd_methods": (2, ';')},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'pdz': input_formats.NetworkInput(
        name = 'PDZBase',
        separator = None,
        id_col_a = 0,
        id_col_b = 2,
        id_type_a = 'uniprot',
        id_type_b = 'uniprot',
        entity_type_a = 'protein',
        entity_type_b = 'protein',
        is_directed = 1,
        sign = False,
        ncbi_tax_id = {
            'col': 7,
            'include': {9606},
        },
        input = 'pdzbase_interactions',
        references = 8,
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}
    ),
    'signor': input_formats.NetworkInput(
        name = 'SIGNOR',
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        # only direct interactions
        positive_filters = [(10, True)],
        # exclude TF-target interactions
        negative_filters = [(7, 'transcriptional regulation')],
        entity_type_a = "protein",
        entity_type_b = "protein",
        ncbi_tax_id = {
            'col': 8,
            'dict': {
                '9606;9606': 9606,
                '9606': 9606
            },
        },
        is_directed = (
            6,
            {
                'up-regulates',
                'up-regulates activity',
                'up-regulates quantity',
                'up-regulates quantity by stabilization',
                'up-regulates quantity by expression',
                'down-regulates',
                'down-regulates activity',
                'down-regulates quantity by destabilization',
                'down-regulates quantity',
                'down-regulates quantity by repression',
            }
        ),
        sign = (
            6,
            {
                'up-regulates',
                'up-regulates activity',
                'up-regulates quantity',
                'up-regulates quantity by stabilization',
                'up-regulates quantity by expression',
            },
            {
                'down-regulates',
                'down-regulates activity',
                'down-regulates quantity by destabilization',
                'down-regulates quantity',
                'down-regulates quantity by repression',
            }
        ),
        input = 'signor.signor_interactions',
        references = (9, ";"),
        header = False,
        extra_edge_attrs = {"signor_mechanism": (7, ';')},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {},
    ),
    'adhesome': input_formats.NetworkInput(
        name = "Adhesome",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "genesymbol",
        id_type_b = "genesymbol",
        entity_type_a = "protein",
        entity_type_b = "protein",
        ncbi_tax_id = 9606,
        is_directed = (2, ('+', '_')),
        sign = (2, '+', '_'),
        input = 'adhesome.adhesome_interactions',
        references = 4,
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {},
    ),
    'icellnet': input_formats.NetworkInput(
        name = "ICELLNET",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        ncbi_tax_id = 9606,
        is_directed = True,
        sign = None,
        input = 'icellnet.icellnet_interactions',
        references = 6,
        resource = 5,
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {},
    ),
}

# synonym
activity_flow = pathway

"""
Pathway (activity flow) resources without literature references.
"""
pathway_noref = {
    'kegg': input_formats.NetworkInput(
        name = "KEGG",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        ncbi_tax_id = 9606,
        is_directed = (2, ('activation', 'inhibition')),
        sign = (2, 'activation', 'inhibition'),
        input = 'kegg.kegg_interactions',
        references = False,
        must_have_references = False,
        header = False,
        positive_filters = [
            (5, True), # is_direct
        ],
        negative_filters = [
            (6, True), # transcriptional
        ],
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {},
    ),
    'wang': input_formats.NetworkInput(
        name = "Wang",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "genesymbol",
        id_type_b = "genesymbol",
        entity_type_a = "protein",
        entity_type_b = "protein",
        ncbi_tax_id = 9606,
        is_directed = (2, ('+', '-')),
        sign = (2, '+', '-'),
        input = 'wang_interactions',
        references = False,
        must_have_references = False,
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {},
    ),
    'kegg-medicus': input_formats.NetworkInput(
        name = "KEGG-MEDICUS",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        ncbi_tax_id = 9606,
        is_directed = (5, ('stimulation', 'inhibition')),
        sign = (5, 'stimulation', 'inhibition'),
        input = 'kegg.kegg_medicus_interactions',
        references = False,
        must_have_references = False,
        header = False,
        positive_filters = [
            (4, 'post_translational'),
        ],
        negative_filters = [
            (5, ['missing', 'enzyme_enzyme']),
        ],
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {},
    ),
}


pathway_all = dict(copy.deepcopy(pathway), **copy.deepcopy(pathway_noref))


pathway_bad = {
    'laudanna_effects': input_formats.NetworkInput(
        name = "Laudanna-effects",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "genesymbol",
        id_type_b = "genesymbol",
        entity_type_a = "protein",
        entity_type_b = "protein",
        ncbi_tax_id = 9606,
        is_directed = (2, ('activation', 'inhibition')),
        sign = (2, 'activation', 'inhibition'),
        input = 'get_laudanna_effects',
        references = False,
        must_have_references = False,
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {},
    ),
    'laudanna_directions': input_formats.NetworkInput(
        name = "Laudanna-directions",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "genesymbol",
        id_type_b = "genesymbol",
        entity_type_a = "protein",
        entity_type_b = "protein",
        ncbi_tax_id = 9606,
        is_directed = True,
        input = 'get_laudanna_directions',
        references = False,
        must_have_references = False,
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {},
    ),
}

'''
Interaction databases included in OmniPath.
These are subsets of the named databases, having
only low throughput, manually curated, undirected
interactions with literature references.
'''
interaction = {
    'biogrid': input_formats.NetworkInput(
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
    'ccmap': input_formats.NetworkInput(
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
    'mppi': input_formats.NetworkInput(
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
        input = 'mppi.mppi_interactions',
        references = (0, ";"),
        ncbi_tax_id = 9606,
        extra_edge_attrs = {"mppi_evidences": (1, ";")},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'dip': input_formats.NetworkInput(
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
    'netpath': input_formats.NetworkInput(
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
        input = 'netpath.netpath_interactions',
        references = (4, ";"),
        ncbi_tax_id = 9606,
        extra_edge_attrs = {
            "netpath_methods": (5, ";"),
            "netpath_type": (6, ";"),
            "netpath_pathways": (7, ';')
        },
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'innatedb': input_formats.NetworkInput(
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
    'alz': input_formats.NetworkInput(
        name = 'AlzPathway',
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = False,
        sign = False,
        input = urls.urls['alzpathway']['url'],
        references = (8, ";"),
        ncbi_tax_id = 9606,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'matrixdb': input_formats.NetworkInput(
        name = "MatrixDB",
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = False,
        sign = False,
        input = 'matrixdb.matrixdb_interactions',
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
    'psite': input_formats.NetworkInput(
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
        input = 'phosphosite.phosphosite_interactions_curated',
        references = (5, ";"),
        ncbi_tax_id = 9606,
        extra_edge_attrs = {"psite_evidences": (4, ";")},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}
    ),
    'depod': input_formats.NetworkInput(
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
        input = 'depod.depod_interactions',
        references = (2, "|"),
        ncbi_tax_id = 9606,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'lmpid': input_formats.NetworkInput(
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
    'phelm': input_formats.NetworkInput(
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
        input = 'phosphoelm.phosphoelm_interactions',
        references = (2, ';'),
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'elm': input_formats.NetworkInput(
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
            'A': {'col': 13, 'include': {9606}},
            'B': {'col': 14, 'include': {9606}},
        },
        input = 'elm_interactions',
        references = 12,
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}
    ),
    'domino': input_formats.NetworkInput(
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
    'dbptm': input_formats.NetworkInput(
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
        input = 'dbptm.dbptm_interactions',
        references = (2, ';'),
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {},
        must_have_references = True
    ),
    'hprd_p': input_formats.NetworkInput(
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
        input = 'hprd.hprd_interactions',
        references = (10, ','),
        header = False,
        extra_edge_attrs = {'hprd_mechanism': 8},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}
    ),
    'protmapper': input_formats.NetworkInput(
        name = "ProtMapper",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = 'uniprot',
        id_type_b = 'uniprot',
        entity_type_a = 'protein',
        entity_type_b = 'protein',
        is_directed = 1,
        sign = False,
        ncbi_tax_id = 9606,
        input = 'protmapper.protmapper_interactions',
        references = 3,
        resource = 2,
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {},
        must_have_references = True,
    ),
    'kea': input_formats.NetworkInput(
        name = 'KEA',
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = 'uniprot',
        id_type_b = 'uniprot',
        entity_type_a = 'protein',
        entity_type_b = 'protein',
        is_directed = 1,
        sign = False,
        ncbi_tax_id = 9606,
        input = 'kea.kea_interactions',
        references = 4,
        resource = 5,
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {},
        must_have_references = True,
    ),
    'iptmnet': input_formats.NetworkInput(
        name = 'iPTMnet',
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = 'uniprot',
        id_type_b = 'uniprot',
        entity_type_a = 'protein',
        entity_type_b = 'protein',
        is_directed = 1,
        sign = False,
        ncbi_tax_id = 9606,
        input = 'iptmnet.iptmnet_interactions',
        references = 8,
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {},
        must_have_references = True,
    ),
}

# synonym
enzyme_substrate = ptm

'''
Other PTM datasets which are not used because the lack of
references.
'''
ptm_misc = {
    'psite_noref': input_formats.NetworkInput(
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
        input = 'phosphosite.phosphosite_interactions_noref',
        references = False,
        extra_edge_attrs = {"psite_evidences": (4, ";")},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}
    ),
    'ppoint': input_formats.NetworkInput(
        name = "PhosphoPoint",
        separator = ";",
        id_col_a = 0,
        id_col_b = 2,
        id_type_a = "genesymbol",
        id_type_b = "genesymbol",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = True,
        header = True,
        ncbi_tax_id = 9606,
        input = 'phosphopoint_interactions',
        references = False,
        sign = False,
        extra_edge_attrs = {"phosphopoint_category": 4},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}
    ),
    'pnetworks': input_formats.NetworkInput(
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
        input = 'phosphonetworks.phosphonetworks_interactions',
        references = False,
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}
    ),
    'mimp': input_formats.NetworkInput(
        name = 'MIMP',
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
        input = 'mimp.mimp_interactions',
        references = False,
        resource = 2,
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}
    ),
    'li2012': input_formats.NetworkInput(
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
        input = 'li2012.li2012_interactions',
        references = False,
        ncbi_tax_id = 9606,
        extra_edge_attrs = {
            'li2012_mechanism': 3,
            'li2012_route': 2,
        },
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}
    ),

}

ptm_misc['protmapper'] = copy.deepcopy(ptm['protmapper'])
ptm_misc['protmapper'].must_have_references = False

# synonym
ptm_noref = ptm_misc

ptm_all = copy.deepcopy(ptm_misc)
ptm_all.update(ptm)

extra_directions = copy.deepcopy(ptm_misc)
extra_directions.update(copy.deepcopy(pathway_noref))
extra_directions['acsn'] = copy.deepcopy(reaction_pc['acsn'])
'''
Interaction databases not included in OmniPath.
These were omitted because lack of references,
or because we could not separate the low throughput,
manually curated interactions.
'''
interaction_misc = {
    'intact': input_formats.NetworkInput(
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
        input = 'intact.intact_interactions',
        references = 4,
        ncbi_tax_id = 9606,
        extra_edge_attrs = {"intact_methods": 5},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {},
    ),
    'hippie': input_formats.NetworkInput(
        name = "HIPPIE",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = False,
        sign = False,
        input = 'hippie.hippie_interactions',
        references = 4,
        ncbi_tax_id = 9606,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'biogrid': input_formats.NetworkInput(
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
    'hi2': input_formats.NetworkInput(
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
        input = 'huri.rolland_hi_ii_14',
        references = False,
        header = False,
        extra_edge_attrs = {'hi2_numof_screens': 4},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'lit13': input_formats.NetworkInput(
        name = "Lit-BM-13",
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
        input = 'huri.lit_bm_13_interactions',
        references = False,
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'lit17': input_formats.NetworkInput(
        name = "Lit-BM-17",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = False,
        sign = False,
        ncbi_tax_id = 9606,
        references = (2, ';'),
        input = 'huri.lit_bm_17_interactions',
        header = False,
        extra_edge_attrs = {
            'mentha_score': 3,
        },
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'huri': input_formats.NetworkInput(
        name = "HuRI",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = False,
        sign = False,
        ncbi_tax_id = 9606,
        input = 'huri.huri_interactions',
        references = False,
        header = True,
        extra_edge_attrs = {
            'huri_score': 4,
        },
        extra_node_attrs_a = {},
        extra_node_attrs_b = {},
    ),
    'hi_union': input_formats.NetworkInput(
        name = "HI-union",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = False,
        sign = False,
        ncbi_tax_id = 9606,
        input = 'huri.hi_union_interactions',
        references = False,
        header = True,
        extra_edge_attrs = {
            'huri_score': 4,
        },
        extra_node_attrs_a = {},
        extra_node_attrs_b = {},
    ),
    'yu2011': input_formats.NetworkInput(
        name = "Yu2011",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = False,
        sign = False,
        ncbi_tax_id = 9606,
        input = 'huri.yu2011_interactions',
        references = False,
        header = True,
        extra_edge_attrs = {
            'huri_score': 4,
        },
        extra_node_attrs_a = {},
        extra_node_attrs_b = {},
    ),
    'yang2016': input_formats.NetworkInput(
        name = "Yang2016",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = False,
        sign = False,
        ncbi_tax_id = 9606,
        input = 'huri.yang2016_interactions',
        references = False,
        header = True,
        extra_edge_attrs = {
            'huri_score': 4,
        },
        extra_node_attrs_a = {},
        extra_node_attrs_b = {},
    ),
    'cpdb': input_formats.NetworkInput(
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
        extra_node_attrs_b = {}),
}

interaction_misc['acsn'] = copy.deepcopy(reaction_pc['acsn'])


interaction_deprecated = {
    'hi3_local': input_formats.NetworkInput(
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
        extra_node_attrs_b = {}
    ),
    'hi3_local_2': input_formats.NetworkInput(
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
        input = 'huri.vidal_hi_iii',
        references = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {},
        must_have_references = False,
        input_args = {'fname': '/home/denes/Documents/pw/data/hi3-2.3.tsv'}
    ),
    'hi3': input_formats.NetworkInput(
        name = "HI-III",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = False,
        sign = False,
        ncbi_tax_id = 9606,
        input = 'huri.hi_iii',
        references = False,
        header = True,
        extra_edge_attrs = {
            'hi3_score': 5,
        },
        extra_node_attrs_a = {},
        extra_node_attrs_b = {},
    ),
    'hsn': input_formats.NetworkInput(
        name = "HumanSignalingNetwork",
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
        extra_node_attrs_b = {}
    ),
}


interaction_htp = {
    'intact': input_formats.NetworkInput(
        name = "IntAct",
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = False,
        sign = False,
        input = 'intact_interactions',
        references = 4,
        ncbi_tax_id = 9606,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {},
        input_args = {'miscore': 0.0}
    ),
    'biogrid': input_formats.NetworkInput(
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
        input_args = {
            'htp_limit': None,
            'ltp': False,
        },
    ),
    'dip': input_formats.NetworkInput(
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
        input_args = {
            'core_only': False,
            'small_scale_only': False,
        }
    ),
    'ccmap': input_formats.NetworkInput(
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
        extra_node_attrs_b = {},
    ),
    'innatedb': input_formats.NetworkInput(
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
        extra_node_attrs_b = {}
    ),
    'matrixdb': input_formats.NetworkInput(
        name = "MatrixDB",
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = False,
        sign = False,
        input = 'matrixdb.matrixdb_interactions',
        references = (2, "|"),
        ncbi_tax_id = 9606,
        extra_edge_attrs = {"matrixdb_methods": (3, '|')},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}
    ),
    'hprd': input_formats.NetworkInput(
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
        input = 'hprd.hprd_interactions_htp',
        references = (7, ','),
        header = False,
        extra_edge_attrs = {'hprd_methods': (6, ';')},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}
    ),
    'mppi': input_formats.NetworkInput(
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
        input = 'mppi.mppi_interactions',
        references = (0, ";"),
        ncbi_tax_id = 9606,
        extra_edge_attrs = {"mppi_evidences": (1, ";")},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}
    ),
}

'''
Transcriptional regulatory interactions.
'''
transcription_onebyone = {
    'abs': input_formats.NetworkInput(
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
        input = 'abs.abs_interactions',
        interaction_type = 'transcriptional',
        references = False,
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {},
    ),
    'encode_dist': input_formats.NetworkInput(
        name = "ENCODE-distal",
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
        interaction_type = 'transcriptional',
        references = False,
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {},
    ),
    'encode_prox': input_formats.NetworkInput(
        name = "ENCODE-proximal",
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
        interaction_type = 'transcriptional',
        references = False,
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {},
    ),
    'pazar': input_formats.NetworkInput(
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
        interaction_type = 'transcriptional',
        references = 2,
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {},
    ),
    'htri': input_formats.NetworkInput(
        name = "HTRIdb",
        separator = None,
        id_col_a = 1,
        id_col_b = 3,
        id_type_a = "genesymbol",
        id_type_b = "genesymbol",
        entity_type_a = "protein",
        entity_type_b = "protein",
        is_directed = True,
        sign = False,
        ncbi_tax_id = 9606,
        input = 'get_htri',
        interaction_type = 'transcriptional',
        references = 4,
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {},
    ),
    'oreganno': input_formats.NetworkInput(
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
        interaction_type = 'transcriptional',
        references = 2,
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {},
    ),
    'signor': input_formats.NetworkInput(
        name = "SIGNOR",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        # only direct TF-target interactions
        positive_filters = [
            (10, True),
            (
                7,
                {
                    'transcriptional regulation',
                    'transcriptional activation',
                    'transcriptional repression',
                }
            ),
        ],
        entity_type_a = "protein",
        entity_type_b = "protein",
        ncbi_tax_id = {'col': 8,
                   'dict': {
                       '9606;9606': 9606,
                       '9606': 9606
                   }},
        is_directed = True,
        sign = (
            6,
            {
                'up-regulates',
                'up-regulates activity',
                'up-regulates quantity',
                'up-regulates quantity by stabilization',
                'up-regulates quantity by expression',
            },
            {
                'down-regulates',
                'down-regulates activity',
                'down-regulates quantity by destabilization',
                'down-regulates quantity',
                'down-regulates quantity by repression',
            }
        ),
        input = 'signor.signor_interactions',
        references = (9, ";"),
        header = True,
        extra_edge_attrs = {"signor_mechanism": (7, ';')},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {},
        interaction_type = 'transcriptional',
    ),
    'kegg': input_formats.NetworkInput(
        name = "KEGG",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        ncbi_tax_id = 9606,
        is_directed = (2, ('repression', 'expression')),
        sign = (2, 'expression', 'repression'),
        input = 'kegg.kegg_interactions',
        references = False,
        must_have_references = False,
        header = False,
        positive_filters = [
            (5, True), # is_direct
            (6, True), # transcriptional
        ],
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {},
        interaction_type = 'transcriptional',
    ),
    'kegg-medicus': input_formats.NetworkInput(
        name = "KEGG-MEDICUS",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "uniprot",
        id_type_b = "uniprot",
        entity_type_a = "protein",
        entity_type_b = "protein",
        ncbi_tax_id = 9606,
        is_directed = (5, ('stimulation', 'inhibition')),
        sign = (5, 'stimulation', 'inhibition'),
        input = 'kegg.kegg_medicus_interactions',
        references = False,
        must_have_references = False,
        header = False,
        positive_filters = [
            (4, 'transcriptional'),
        ],
        negative_filters = [
            (5, ['missing', 'enzyme_enzyme']),
        ],
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {},
        interaction_type = 'transcriptional',
    ),
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
transcription_dorothea = {
    'dorothea': input_formats.NetworkInput(
        name = "DoRothEA",
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
        input = 'dorothea.get_dorothea',
        interaction_type = 'transcriptional',
        resource = (12, ','),
        references = (13, ','),
        header = False,
        extra_edge_attrs = {
            'dorothea_curated': 4,
            'dorothea_chipseq': 5,
            'dorothea_tfbs':    6,
            'dorothea_coexp':   7,
            'dorothea_level':   3,
            'dorothea_kegg_pathways': (14, '|'),
        },
        extra_node_attrs_a = {},
        extra_node_attrs_b = {},
        must_have_references = False,
    ),
}

# synonyms
dorothea = transcription_dorothea
tfregulons = transcription_dorothea

# all transcriptional regulation resources
transcription = {}
transcription.update(copy.deepcopy(transcription_onebyone))
transcription.update(copy.deepcopy(transcription_dorothea))

'''
Old transctiptional regulation input formats.
Should not be used.
'''
transcription_deprecated = {
    'oreganno_old': input_formats.NetworkInput(
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
        interaction_type = 'transcriptional',
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
    'mir2dis': input_formats.NetworkInput(
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
        interaction_type = 'post_transcriptional',
        references = None,
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'mirdeath': input_formats.NetworkInput(
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
        interaction_type = 'post_transcriptional',
        references = 3,
        header = False,
        extra_edge_attrs = {'mirdeathdb_function': 4},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'ncrdeath': input_formats.NetworkInput(
        name = "ncRDeathDB",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "mirbase",
        id_type_b = "genesymbol",
        entity_type_a = "mirna",
        entity_type_b = "protein",
        is_directed = True,
        sign = False,
        ncbi_tax_id = {
            'col': 6,
            'include': {9606},
        },
        input = 'ncrdeathdb_interactions',
        interaction_type = 'post_transcriptional',
        references = 5,
        header = False,
        positive_filters = [(2, 'miRNA')],
        negative_filters = [
            (5, 'prediction'),
            (1, {None}),
            (0, {None}),
        ],
        extra_edge_attrs = {
            'ncrdeathdb_pathway': 3,
            'ncrdeathdb_effect': 4,
        },
        extra_node_attrs_a = {},
        extra_node_attrs_b = {},
    ),
    'mirecords': input_formats.NetworkInput(
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
                'dict': common.swap_dict(taxonomy.phosphoelm_taxids),
                'include': set([9606])
            },
            'B': {
                'col': 4,
                'dict': common.swap_dict(taxonomy.phosphoelm_taxids),
                'include': set([9606])
            }},
        input = 'mirecords_interactions',
        interaction_type = 'post_transcriptional',
        references = 5,
        header = False,
        extra_edge_attrs = {'mirdeathdb_function': 4},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'mirtarbase': input_formats.NetworkInput(
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
                'dict': common.swap_dict(taxonomy.phosphoelm_taxids),
                'include': set([9606])
            },
            'B': {
                'col': 5,
                'dict': common.swap_dict(taxonomy.phosphoelm_taxids),
                'include': set([9606])
            }},
        positive_filters = [(7, 'Functional MTI')],
        input = 'mirtarbase_interactions',
        interaction_type = 'post_transcriptional',
        references = 8,
        header = False,
        extra_edge_attrs = {'mirtarbase_evidence': (6, '//')},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {})
}

tf_mirna = {
    'transmir': input_formats.NetworkInput(
        name = "TransmiR",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "genesymbol",
        id_type_b = "mir-mat-name",
        entity_type_a = "protein",
        entity_type_b = "mirna",
        is_directed = True,
        sign = (2, 'Activation', 'Repression'),
        ncbi_tax_id = 9606,
        input = 'transmir_interactions',
        interaction_type = 'mirna_transcriptional',
        references = 3,
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}
    ),
    'encode': input_formats.NetworkInput(
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
        interaction_type = 'mirna_transcriptional',
        references = None,
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {})
}

lncrna_target = {
    'lncdisease': input_formats.NetworkInput(
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
                'dict': common.swap_dict(taxonomy.taxids),
                'include': set([9606])
            },
        positive_filters = [(2, 'RNA'), (3, 'Protein')],
        input = 'lncdisease_interactions',
        interaction_type = 'lncrna_post_transcriptional',
        references = 6,
        header = False,
        extra_edge_attrs = {'lncrnadisease_mechanism': 4},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}),
    'lncrnadb': input_formats.NetworkInput(
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
                'dict': common.swap_dict(taxonomy.phosphoelm_taxids),
                'include': set([9606])
            },
        positive_filters = [(2, 'protein')],
        input = 'lncrnadb_interactions',
        interaction_type = 'lncrna_post_transcriptional',
        references = 4,
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {}
    ),
    'ncrdeath': input_formats.NetworkInput(
        name = "ncRDeathDB",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "lncrna-genesymbol",
        id_type_b = "genesymbol",
        entity_type_a = "lncrna",
        entity_type_b = "protein",
        is_directed = True,
        sign = False,
        ncbi_tax_id = {
            'col': 6,
            'include': {9606},
        },
        input = 'ncrdeathdb_interactions',
        interaction_type = 'lncrna_post_transcriptional',
        references = 5,
        header = False,
        positive_filters = [(2, 'lncRNA')],
        negative_filters = [
            (5, 'prediction'),
            (1, {None}),
            (0, {None}),
        ],
        extra_edge_attrs = {
            'ncrdeathdb_pathway': 3,
            'ncrdeathdb_effect': 4,
        },
        extra_node_attrs_a = {},
        extra_node_attrs_b = {},
    ),
}

ligand_receptor = {
    'ramilowski2015': input_formats.NetworkInput(
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
        input = 'ramilowski2015.ramilowski_interactions',
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
        data_model = 'ligand_receptor',
    ),
    'kirouac2010': input_formats.NetworkInput(
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
        input = 'kirouac2010.kirouac2010_interactions',
        references = False,
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {},
        mark_source = 'kirouac_ligand',
        mark_target = 'kirouac_receptor',
        data_model = 'ligand_receptor',
    ),
    'hpmr': input_formats.NetworkInput(
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
        data_model = 'ligand_receptor',
    ),
    'cellphonedb': input_formats.NetworkInput(
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
        input = 'cellphonedb.cellphonedb_interactions',
        references = (3, ';'),
        resource = (2, ';'),
        must_have_references = False,
        header = False,
        extra_edge_attrs = {'cellphonedb_type': 4},
        extra_node_attrs_a = {'cellphonedb_type': 5},
        extra_node_attrs_b = {'cellphonedb_type': 6},
        positive_filters = [],
        data_model = 'ligand_receptor',
    ),
    'lrdb': input_formats.NetworkInput(
        name = 'LRdb',
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
        input = 'lrdb.lrdb_interactions',
        references = 3,
        resource = 2,
        must_have_references = False,
        header = False,
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {},
        data_model = 'ligand_receptor',
    ),
    'baccin2019': input_formats.NetworkInput(
        name = 'Baccin2019',
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
        input = 'baccin2019.baccin2019_interactions',
        references = 6,
        resource = 5,
        must_have_references = False,
        header = False,
        negative_filters = [
            (2, 'Incorrect'),
        ],
        extra_edge_attrs = {},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {
            'baccin_category': 4,
            'baccin_location': 3,
        },
        data_model = 'ligand_receptor',
    ),
    'embrace': input_formats.NetworkInput(
        name = 'EMBRACE',
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
        input = 'embrace.embrace_interactions',
        must_have_references = False,
        header = False,
        data_model = 'ligand_receptor',
    ),
    'italk': input_formats.NetworkInput(
        name = 'iTALK',
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
        input = 'italk.italk_interactions',
        must_have_references = False,
        header = False,
        data_model = 'ligand_receptor',
    ),
}

small_molecule_protein = {
    'signor': input_formats.NetworkInput(
        name = "SIGNOR",
        separator = None,
        id_col_a = 0,
        id_col_b = 1,
        id_type_a = "pubchem-cid",
        id_type_b = "uniprot",
        # only direct interactions
        positive_filters = [(10, True), (6, 'chemical')],
        # exclude TF-target interactions
        negative_filters = [(7, 'transcriptional regulation')],
        entity_type_a = "small_molecule",
        entity_type_b = "protein",
        ncbi_tax_id = {'col': 8,
                   'dict': {
                       '9606;9606': 9606,
                       '9606': 9606
                   }},
        is_directed = (
            6,
            [
                'up-regulates',
                'up-regulates activity',
                'up-regulates quantity by stabilization',
                'down-regulates',
                'down-regulates activity',
                'down-regulates quantity by destabilization',
            ]
        ),
        sign = (
            6,
            [
                'up-regulates',
                'up-regulates activity',
                'up-regulates quantity by stabilization',
            ],
            [
                'down-regulates',
                'down-regulates activity',
                'down-regulates quantity by destabilization',
            ]
        ),
        input = 'signor.signor_interactions',
        references = (9, ";"),
        header = True,
        extra_edge_attrs = {"signor_mechanism": (7, ';')},
        extra_node_attrs_a = {},
        extra_node_attrs_b = {},
    ),
}

ligand_receptor['guide2pharma'] = copy.deepcopy(pathway['guide2pharma'])
ligand_receptor['guide2pharma'].data_model = 'ligand_receptor'
ligand_receptor['icellnet'] = copy.deepcopy(pathway['icellnet'])
ligand_receptor['icellnet'].must_have_references = False
ligand_receptor['icellnet'].data_model = 'ligand_receptor'

pathway['hpmr'] = copy.deepcopy(ligand_receptor['hpmr'])
pathway['hpmr'].data_model = 'activity_flow'
pathway['hpmr'].must_have_references = True
pathway['hpmr'].positive_filters = []
pathway['cellphonedb'] = copy.deepcopy(ligand_receptor['cellphonedb'])
pathway['cellphonedb'].must_have_references = True
pathway['cellphonedb'].data_model = 'activity_flow'
pathway['ramilowski2015'] = copy.deepcopy(ligand_receptor['ramilowski2015'])
pathway['ramilowski2015'].must_have_references = True
pathway['ramilowski2015'].data_model = 'activity_flow'
pathway['lrdb'] = copy.deepcopy(ligand_receptor['lrdb'])
pathway['lrdb'].data_model = 'activity_flow'
pathway['lrdb'].must_have_references = True
pathway['baccin2019'] = copy.deepcopy(ligand_receptor['baccin2019'])
pathway['baccin2019'].data_model = 'activity_flow'
pathway['baccin2019'].must_have_references = True

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
    'negatome': input_formats.NetworkInput(
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

biocarta = input_formats.NetworkInput(
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

nci_pid = input_formats.NetworkInput(
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

reactome = input_formats.NetworkInput(
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

gdsc_comp_target = input_formats.NetworkInput(
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
