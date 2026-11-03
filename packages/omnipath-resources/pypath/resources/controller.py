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

from __future__ import annotations

"""
Highest level resource management API.
"""

from future.utils import iteritems

from typing import Iterable, Literal

import json
import os
import copy
import importlib as imp
import itertools
import functools

import pandas as pd

import pypath.share.session as session_mod
import pypath.share.common as common
import pypath.internals.resource as resource_base
import pypath.resources._network as netres
from . import licenses as licenses


class ResourceController(session_mod.Logger):
    """
    Resource controller is aimed to be the central part of pypath
    communication with resources.

    14.01.2020: the initial step for resource controller development:
        used for /info page generating for the server.
    14.02.2020: storing and reading enzyme-substrate resource definitions
        from the JSON; class inherits from session.Logger
    """

    def __init__(
            self,
            resource_info_path = (
                session_mod.session().module_root,
                'resources',
                'data',
                'resources.json',
            ),
            use_package_path = False,
        ):

        session_mod.Logger.__init__(self, name = 'resource_controller')

        self.data = None

        if use_package_path:

            resource_info_path = (
                (
                    os.path.dirname(os.path.abspath(__file__)),
                ) +
                resource_info_path
            )

        self.resource_info_path = os.path.join(*resource_info_path)

        self._log(
            'Loading resource information from '
            'JSON file: %s' % self.resource_info_path
        )

        self.update()


    def export_licenses(self, path='licenses.tsv'):
        """
        Exports all resources' license information as a TSV file
        """

        df = pd.DataFrame(
            [
                [
                    k,
                    v.name,
                    v.full_name,
                    v.purpose.level,
                    v.attrib.level,
                    v.sharing.level,
                    v.url,
                ]
                for k, v in self.licenses.items()
            ],
            columns=[
                'resource',
                'name',
                'full_name',
                'purpose',
                'attrib',
                'sharing',
                'url',
            ]
        )

        df.to_csv(path, sep='\t', index=False)


    def reload(self):

        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)


    def update(self, path = None, force = False, remove_old = False):
        """
        Reads resource information from a JSON file.

        :arg str,NoneType path:
            Path to a JSON file with resource information. By default the
            path in py:attr:``resource_info_path`` used which points by
            default to the built in resource information file.
        :arg bool force:
            Read the file again even if no new path provided and it has been
            read already.
        :arg bool remove_old:
            Remove old data before reading. By default the data will be
            updated with the contents of the new file potentially overwriting
            identical keys in the old data.
        """

        if self.data and not path and not force:

            return

        if not self.data or remove_old:

            self.data = {}

        path = path or self.resource_info_path

        try:

            with open(path, encoding = 'utf-8') as json_file:

                resources_data = json.load(json_file)
                self.data = resources_data
                self._log(
                    'Resource information has been read from `%s`.' % path
                )
                self.update_licenses()

        except IOError:

            self._console(
                'File %s with resources information cannot be accessed. '
                'Check the name of the file.' % path
            )


    def update_licenses(self):

        self.license_db = licenses.Licenses()
        self.licenses = {}
        self.synonyms = {}
        self.secondary = {}

        self._log('Updating resource license information.')

        for res, res_data in iteritems(self.data):

            if 'license' in res_data:

                res_data['license'] = self.license_db[res_data]
                self.licenses[res] = res_data['license']

                for synonym in res_data.get('synonyms', ()):

                    self.licenses[synonym] = res_data['license']
                    self.synonyms[synonym] = res

                if 'components' in res_data:

                    self.secondary[res] = set(res_data['components'])

            else:

                self._log(f'No license for resource `{res}`.')


    def __getitem__(self, key):

        return self.resource(key)


    def resource(self, name):

        return self._get(name, dct = self.data)


    def name(self, name):

        if name in self.synonyms:

            name = self.synonyms[name]

        return name


    @functools.cache
    def secondary_resources(self, name, postfix = False):
        """
        Args:
            name:
                Name of a composite resource.
            postfix:
                Append the name of the primary resource to the secondary,
                separated by an underscore, e.g. "TFactS_CollecTRI".
        """

        name = self.name(name)

        secondary = self.secondary.get(name, set())

        if postfix:

            secondary = {f'{sec}_{name}' for sec in secondary}

        return secondary


    def _get(self, name, dct):

        if name in dct:

            return dct[name]

        elif name in self.synonyms:

            name = self.synonyms[name]
            return dct[name]

        elif '_' in name:

            name = name.split('_', maxsplit = 1)[0]
            return self._get(name, dct)

        else:

            self._log('Could not find resource `%s`.' % name)


    def license(self, name):

        return self._get(name, dct = self.licenses)


    def license_filter(
            self,
            resources: list | dict,
            purpose: Literal[
                'academic',
                'commercial',
                'for-profit',
                'non-profit',
                'ignore',
            ] | None = None,
            sharing: Literal[
                'alike',
                'free',
                'noderiv',
                'noshare',
                'share',
                'deriv',
                'ignore',
            ] | None = None,
            attrib: Literal[
                'attrib',
                'free',
                'noattrib',
                'composite',
                'ignore',
            ] | None = None,
        ) -> list | dict:
        """
        Filters a list of resources by their license.
        """

        self.add_resource_attrs(resources)

        return common.compr(
            obj = resources,
            filter = lambda r: r.license.enables(purpose, sharing, attrib),
        )


    def add_resource_attrs(
            self,
            resources: dict | Iterable[resource_base.AbstractResource],
        ) -> None:
        """
        Adds resource attributes to a list of resources.

        It modifies the instances in-place, returns nothing.
        """

        _ = common.compr(
            resources,
            lambda r: setattr(r, 'resource_attrs', self.resource(r.name)),
        )


    def collect(self, data_type):

        resource_cls = getattr(
            resource_base,
            '%sResource' % (
                ''.join(n.capitalize() for n in data_type.split('_'))
            )
        )

        result = []

        for name, attrs in iteritems(self.data):

            if 'inputs' in attrs and data_type in attrs['inputs']:

                args = copy.deepcopy(attrs['inputs'][data_type])
                args['resource_attrs'] = attrs
                if 'name' not in args:
                    args['name'] = name

                result.append(
                    resource_cls(**args)
                )

        return result


    def collect_enzyme_substrate(self):

        return self.collect('enzyme_substrate')


    def collect_network(
            self,
            datasets: Iterable[
                Literal[
                    'pathway',
                    'pathway_noref',
                    'pathway_all',
                    'activity_flow',
                    'mirna_target',
                    'dorothea',
                    'tfregulons',
                    'omnipath',
                    'reaction_pc',
                    'enzyme_substrate',
                    'extra_directions',
                    'small_molecule_protein',
                    'tf_mirna',
                    'pathwaycommons',
                    'pathwaycommons_transcription',
                    'interaction',
                    'interaction_htp',
                    'interaction_misc',
                    'ligand_receptor',
                    'lncrna_target',
                    'transcription_onebyone',
                    'transcription_dorothea',
                    'ptm',
                    'ptm_noref',
                    'ptm_all',
                    'reaction',
                    'reaction_misc',
                    'negative',
                ],
            ] | None = 'pathway',
            interaction_types: Iterable[
                Literal[
                    'post_translational',
                    'transcriptional',
                    'small_molecule_protein',
                    'post_transcriptional',
                ],
            ] | None = 'post_translational',
            data_models: Iterable[
                Literal[
                    'activity_flow',
                    'interaction',
                    'enzyme_substrate',
                    'process_description',
                    'ligand_receptor',
                    'drug_target',
                ],
            ] | None = 'activity_flow',
            license_purpose: Literal[
                'academic',
                'commercial',
                'for-profit',
                'non-profit',
                'ignore',
            ] = 'ignore',
            license_sharing: Literal[
                'alike',
                'free',
                'noderiv',
                'noshare',
                'share',
                'deriv',
                'ignore',
            ] = 'ignore',
            license_attrib: Literal[
                'attrib',
                'free',
                'noattrib',
                'composite',
                'ignore',
            ] = 'ignore',
            **kwargs
        ) -> dict:
        """
        Collect network (interaction) resource definitions.

        Args:
            interaction_types:
                Include only these interaction types.
            data_models:
                Inclde only these data models.
            datasets:
                Process only these datasets. Note: there are many synonyms
                and overlaps among datasets. In addition, the overlaps might
                apply slightly different settings for the same resource, e.g.
                in `pathway`, interactions must have literature references,
                while in `pathway_noref` the same resources might allow
                interactions without references. The safest is to process only
                one dataset at a time and load them into the `Network` object
                sequentially.
            license_purpose:
                Do not include the resources that are not legally compatible
                with the defined purpose.
            license_sharing:
                Include only resources that allow the desired redistribution
                conditions. E.g. "deriv" means that the resources must allow
                the sharing of their derivative (altered) versions.
            license_attrib:
                Include only resources that allow the desired level of
                attribution. E.g. "noattrib" means that you can use the
                resource without even mentioning who created it.
            kwargs:
                Custom filters. Names should be attributes of the resource
                or the `NetworkInput` object. The special key `__resource__`
                can be used to refer to the whole `NetworkResource` object.
                For simple values, the test is equality, for arrays incidence,
                while custom callables can be provided for more flexible
                filters.
        """

        interaction_types = common.to_set(interaction_types)
        data_models = common.to_set(data_models)
        datasets = common.to_set(datasets)

        kwargs = {
            k: v if callable(v) else lambda x: x in common.to_set(v)
            for k, v in kwargs.items()
        }

        resources = itertools.chain(*(
            getattr(netres, dset).items()
            for dset in datasets
        ))

        resources = {
            key: res
            for key, res in resources
            if (
                not interaction_types or
                res.interaction_type in interaction_types
            ) and
            (
                not datasets or
                res.data_model in data_models
            ) and
            all(
                fltr(
                    res
                        if key == '__resource__' else
                    getattr(res, getattr(res.networkinput, key))
                )
                for key, fltr in kwargs.items()
            )
        }

        resources = self.license_filter(
            resources,
            purpose = license_purpose,
            sharing = license_sharing,
            attrib = license_attrib,
        )

        return resources

    # synonym
    collect_interaction = collect_network
