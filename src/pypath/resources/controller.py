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

import json
import os

import pypath.share.session as session_mod
import pypath.share.common as common
import pypath.internals.resource as resource_base


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
                'resources_descriptions.json',
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

        except IOError:

            self._console(
                'File %s with resources information cannot be accessed. '
                'Check the name of the file.' % path
            )


    def collect_resource_definitions(self, data_type):

        


    def collect_enzyme_substrate_resource_definitions(self):

        
