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

import copy

import pypath.internals.resource as resource
import pypath.resources.data_formats as data_formats
import pypath.share.session as session_mod
import pypath.share.settings as settings

_logger = session_mod.Logger(name = 'network_resources')
_log = _logger._log

_data_models = {
    'interaction': 'interaction',
    'interaction_misc': 'interaction',
    'interaction_htp': 'interaction',
    'ligand_receptor': 'ligand_receptor',
    'ligrecextra': 'ligand_receptor',
    'pathway': 'activity_flow',
    'pathway_all': 'activity_flow',
    'pathway_noref': 'activity_flow',
    'pathwaycommons': 'activity_flow',
    'activity_flow': 'activity_flow',
    'pathwayextra': 'activity_flow',
    'dorothea': 'activity_flow',
    'collectri': 'activity_flow',
    'tf_target': 'activity_flow',
    'transcription': 'activity_flow',
    'transcription_dorothea': 'activity_flow',
    'transcription_collectri': 'activity_flow',
    'transcription_onebyone': 'activity_flow',
    'pathwaycommons_transcription': 'activity_flow',
    'tfregulons': 'activity_flow',
    'mirna_target': 'activity_flow',
    'lncrna_target': 'activity_flow',
    'tf_mirna': 'activity_flow',
    'enzyme_substrate': 'enzyme_substrate',
    'ptm': 'enzyme_substrate',
    'ptm_all': 'enzyme_substrate',
    'ptm_misc': 'enzyme_substrate',
    'ptm_noref': 'enzyme_substrate',
    'kinaseextra': 'enzyme_substrate',
    'reaction': 'process_description',
    'reaction_misc': 'process_description',
    'reaction_pc': 'process_description',
}


DATASET_PRIORITIES = {
    'omnipath': 0,
    'directionextra': 100,
}



def dataset_priority(dataset: str, default: int = 50) -> int:

    return DATASET_PRIORITIES.get(dataset, default)


def choose_dataset(dataset_a: str, dataset_b: str) -> str:

    priority_a = dataset_priority(dataset_a)
    priority_b = dataset_priority(dataset_b)

    return dataset_a if priority_a < priority_b else dataset_b


def _networkinput_to_networkresource(networkinput, data_model = None):

    return resource.NetworkResource(
        name = networkinput.name,
        interaction_type = networkinput.interaction_type,
        networkinput = networkinput,
        data_model = data_model,
    )


def dorothea_expand_levels(resources = None, levels = None):
    """
    In a dictionary of resource definitions creates a separate
    ``NetworkResource`` object for each confidence levels of DoRothEA
    just like each level was a different resource.

    No matter ``resources`` is a ``NetworkResource`` or a dict of network
    resources, returns always a dict of network resources.
    """

    resources = resources or transcription
    levels = levels or settings.get('tfregulons_levels')
    dorothea = {}

    dorothea_original = (
        resources
            if hasattr(resources, 'networkinput') else
        resources['dorothea']
            if 'dorothea' in resources else
        transcription['dorothea']
    )

    for level in levels:

        level_key = 'dorothea_%s' % level

        dorothea[level_key] = copy.deepcopy(dorothea_original)
        dorothea[level_key].name = 'DoRothEA_%s' % level
        dorothea[level_key].networkinput.name = 'DoRothEA_%s' % level
        dorothea[level_key].networkinput.input_args = {'levels': {level}}

    if resources:

        resources = copy.deepcopy(resources)
        _ = resources.pop('dorothea', None)
        resources.update(dorothea)

        return resources

    return dorothea


for dataset_label in dir(data_formats):

    dataset = getattr(data_formats, dataset_label)

    if not isinstance(dataset, dict):

        continue

    new_dataset = resource.NetworkDataset(name = dataset_label)

    for resource_label, input_def in iteritems(dataset):

        if not isinstance(input_def, data_formats.input_formats.NetworkInput):

            continue

        data_model = (
            input_def.data_model or
            (
                _data_models[dataset_label]
                    if dataset_label in _data_models else
                'unknown'
            )
        )

        if (
            data_model == 'unknown' and
            dataset_label not in {'omnipath', 'extra_directions'}
        ):

            _log(
                'Could not find data model for '
                'resource `%s` in set `%s`.' % (
                    input_def.name,
                    dataset_label,
                )
            )

        new_dataset[resource_label] = _networkinput_to_networkresource(
            networkinput = input_def,
            data_model = data_model,
        )

    if new_dataset:

        globals()[dataset_label] = new_dataset


# these we need to re-create to have the data models set correctly
extra_directions = copy.deepcopy(ptm_misc)
extra_directions.update(copy.deepcopy(pathway_noref))
extra_directions['acsn'] = copy.deepcopy(reaction_pc['acsn'])
extra_directions['acsn'].data_model = 'activity_flow'
omnipath = resource.NetworkDataset(name = 'omnipath')
omnipath.update(pathway)
omnipath.update(ptm)
omnipath.update(interaction)
omnipath['hprd'] = interaction_htp['hprd']
