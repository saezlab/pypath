#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#  (Planned for) centrally handling cache for all databases/resources.
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

from future.utils import iteritems

import copy

import pypath.internals.resource as resource
import pypath.resources.data_formats as data_formats
import pypath.share.session as session_mod

_logger = session_mod.Logger(name = 'network_resources')
_log = _logger._log

_data_models = {
    'interaction': 'interaction',
    'interaction_misc': 'interaction',
    'interaction_htp': 'interaction',
    'ligand_receptor': 'ligand_receptor',
    'pathway': 'activity_flow',
    'pathway_all': 'activity_flow',
    'pathway_noref': 'activity_flow',
    'activity_flow': 'activity_flow',
    'dorothea': 'activity_flow',
    'transcription': 'activity_flow',
    'transcription_dorothea': 'activity_flow',
    'transcription_onebyone': 'activity_flow',
    'tfregulons': 'activity_flow',
    'mirna_target': 'activity_flow',
    'lncrna_target': 'activity_flow',
    'tf_mirna': 'activity_flow',
    'enzyme_substrate': 'enzyme_substrate',
    'ptm': 'enzyme_substrate',
    'ptm_all': 'enzyme_substrate',
    'ptm_misc': 'enzyme_substrate',
    'ptm_noref': 'enzyme_substrate',
    'reaction': 'process_description',
    'reaction_misc': 'process_description',
    'reaction_pc': 'process_description',
}


def _networkinput_to_networkresource(networkinput, data_model = None):
    
    return resource.NetworkResource(
        name = networkinput.name,
        interaction_type = networkinput.interaction_type,
        networkinput = networkinput,
        data_model = data_model,
    )


for resource_set_label in dir(data_formats):
    
    resource_set = getattr(data_formats, resource_set_label)
    
    if not isinstance(resource_set, dict):
        
        continue
    
    new_resource_set = {}
    
    for resource_label, input_def in iteritems(resource_set):
        
        if not isinstance(input_def, data_formats.input_formats.NetworkInput):
            
            continue
        
        data_model = (
            input_def.data_model or
            (
                _data_models[resource_set_label]
                    if resource_set_label in _data_models else
                'unknown'
            )
        )
        
        if (
            data_model == 'unknown' and
            resource_set_label not in {'omnipath', 'extra_directions'}
        ):
            
            _log(
                'Could not find data model for '
                'resource `%s` in set `%s`.' % (
                    input_def.name,
                    resource_set_label,
                )
            )
        
        new_resource_set[resource_label] = _networkinput_to_networkresource(
            networkinput = input_def,
            data_model = data_model,
        )
    
    if new_resource_set:
        
        globals()[resource_set_label] = new_resource_set


# these we need to re-create to have the data models set correctly
extra_directions = copy.deepcopy(ptm_misc)
extra_directions.update(copy.deepcopy(pathway_noref))
extra_directions['acsn'] = copy.deepcopy(reaction_pc['acsn'])
extra_directions['acsn'].data_model = 'activity_flow'
omnipath = {}
omnipath.update(pathway)
omnipath.update(ptm)
omnipath.update(interaction)
