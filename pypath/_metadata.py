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

"""
Package metadata (version, authors, etc).
"""

__all__ = ['get_metadata']

import os
import pathlib
import importlib.metadata

import toml
import itertools

_VERSION = '0.16.20'


def get_metadata():
    """
    Basic package metadata.

    Retrieves package metadata from the current project directory or from
    the installed package.
    """

    here = pathlib.Path(__file__).parent
    pyproj_toml = 'pyproject.toml'
    meta = {}

    for project_dir in (here, here.parent, here.parent.parent):

        toml_path = str(project_dir.joinpath(pyproj_toml).absolute())

        if os.path.exists(toml_path):

            pyproject = toml.load(toml_path)

            project = pyproject.get('project')
            project = project or pyproject.get('tool', {}).get('poetry', {})

            packages = set(
                itertools.chain(*(
                    package.values()
                    for package in project.get('packages', [])
                ))
            )

            name = __name__.split('.')[0]
            if name == project.get('name') or name in packages:

                meta = {
                    'name': project.get('name'),
                    'version': project.get('version'),
                    'author': project.get('authors'),
                    'license': project.get('license'),
                    'full_metadata': pyproject,
                }

                break

    if not meta:

        try:

            meta = {
                k.lower(): v for k, v in
                importlib.metadata.metadata(here.name).items()
            }

        except importlib.metadata.PackageNotFoundError:

            pass

    meta['version'] = meta.get('version', None) or _VERSION

    return meta


metadata = get_metadata()
__version__ = metadata.get('version', None)
__author__ = metadata.get('author', None)
__license__ = metadata.get('license', None)
