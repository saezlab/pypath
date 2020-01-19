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

import pypath.session_mod as session_mod
_logger = session_mod.Logger(name='resource_controller')


class ResourcesController:
    """
    Resource controller is aimed to be the central part of pypath communication with resources.
    14.01.2020: the initial step for resource controller development: used for /info page generating for the server.
    """

    def __init__(self, list_resources_path="resources_jsons/resources_descriptions.json", use_package_path=False):
        if use_package_path:
            abs_path = os.path.dirname(os.path.abspath(__file__))
            self.list_resources_path = abs_path + "/" + list_resources_path
        else:
            self.list_resources_path = list_resources_path

        _logger._log("Resources list reference to this json: %s" % self.list_resources_path)

    def get_info_all_resources(self):
        """
        :return: list of of available resources in pypath
        """
        resources_data = []
        try:
            with open(self.list_resources_path) as json_file:
                resources_data = json.load(json_file)
        except IOError:
            _logger._console("File %s with resources information cannot be accessed. Check the name of the file." %
                                str(self.list_resources_path))
        return resources_data
