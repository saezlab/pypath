#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2022
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: Dénes Türei (turei.denes@gmail.com)
#           Nicolàs Palacio
#           Sebastian Lobentanzer
#           Erva Ulusoy
#           Olga Ivanova
#           Ahmet Rifaioglu
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

"""
``pypath`` is a module primarily for building molecular interaction networks
but also with several submodules for accessing, preprocessing and serving
data from various resources.
"""

import sys
import os
import importlib

import pypath._version as _version_mod
import pypath.share.session as _session_mod

_session_mod.new_session()
__version__ = _version_mod.__version__
session = _session_mod.get_session()


def log():
    """
    Browse the current pypath logfile.
    """

    _session_mod.get_log().browse()


class pypath(object):


    __version__ = __version__
    __author__ = _version_mod.__author__

    _disclaimer_text = (
        '\n\t=== d i s c l a i m e r ===\n\n'
        '\tAll data accessed through this module,\n'
        '\teither as redistributed copy or downloaded using the\n'
        '\tprogrammatic interfaces included in the present module,\n'
        '\tare free to use at least for academic research or\n'
        '\teducation purposes.\n'
        '\tPlease be aware of the licenses of all the datasets\n'
        '\tyou use in your analysis, and please give appropriate\n'
        '\tcredits for the original sources when you publish your\n'
        '\tresults. To find out more about data sources please\n'
        '\tlook at `pypath/resources/data/resources.json` or\n'
        '\thttps://omnipathdb.org/info and \n'
        '\t`pypath.resources.urls.urls`.\n\n'
    )


    @classmethod
    def _disclaimer(cls):

        sys.stdout.write(cls._disclaimer_text)
        sys.stdout.flush()


def info():

    pypath._disclaimer()

    _session_mod.get_log().msg(
        (
            '\n'
            '\t- session ID: `%s`\n'
            '\t- working directory: `%s`\n'
            '\t- logfile: `%s`\n'
            '\t- pypath version: %s' % (
                _session_mod.get_session().label,
                os.getcwd(),
                _session_mod.get_log().fname,
                __version__
            )
        ),
        label = 'pypath',
        level = -9,
        wrap = False,
    )
