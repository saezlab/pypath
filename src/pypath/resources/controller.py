#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
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

import json
import os
import copy
import importlib as imp

import pypath.share.session as session_mod
import pypath.share.common as common
import pypath.internals.resource as resource_base
import pypath.resources.licenses as licenses


_logger = session_mod.Logger(name = 'resources.controller')


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
                common.ROOT,
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

            with open(path) as json_file:

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

                if isinstance(res_data['license'], common.basestring):

                    self._update_license(res_data)
                    self.licenses[res] = res_data['license']

                    if 'synonyms' in res_data:

                        for synonym in res_data['synonyms']:

                            self.licenses[synonym] = res_data['license']
                            self.synonyms[synonym] = res

                    if 'components' in res_data:

                        self.secondary[res] = set(res_data['components'])

            else:

                self._log('No license for resource `%s`.' % res)


    def _update_license(self, resource_data):

        license_key = resource_data['license']
        resource_data['license'] = self.license_db[license_key]


    def __getitem__(self, key):

        return self.resource(key)


    def resource(self, name):

        return self._get(name, dct = self.data)


    def name(self, name):

        if name in self.synonyms:

            name = self.synonyms[name]

        return name


    def secondary_resources(self, name):

        name = self.name(name)

        return self.secondary[name] if name in self.secondary else set()


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
