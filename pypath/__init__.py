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
``pypath`` is a Python module building molecular biology databases. With
``pypath`` you can access, process and combine more than 170 public
resources, build molecular interaction networks, include annotations
about the function, localization, structure and other properties of
genes and proteins. ``pypath`` also provides a number of utilities for
ID translation, orthology translation, handling of taxonomy and sequences,
serving the combined databases by a web server, etc.
"""

import sys
import os
import platform

from pypath._metadata import __version__, __author__, __license__
from pypath._metadata import metadata as _metadata
from pypath.share import session as _session

session = _session.session()


def log():
    """
    Browse the current pypath logfile.
    """

    session._logger.browse()


def disclaimer():
    """
    Prints a disclaimer about copyrights of database data.

    Also points to further information about licenses.
    """

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

    sys.stdout.write(_disclaimer_text)
    sys.stdout.flush()


def info(loglevel = -9):
    """
    Prints basic information about the current session.
    """

    session._logger.msg(
        (
            '\n'
            '\t- session ID: `%s`\n'
            '\t- working directory: `%s`\n'
            '\t- logfile: `%s`\n'
            '\t- config: \n\t\t- %s\n'
            '\t- pypath version: %s\n'
            '\t- imported from: `%s`\n'
            '\t- Python version: %s\n'
            '\t- Platform: %s'% (
                session.label,
                os.getcwd(),
                session._logger.fname,
                '\n\t\t- '.join(
                    map(str, session.config._parsed)
                ),
                __version__,
                os.path.dirname(__file__),
                sys.version,
                platform.platform(),
            )
        ),
        label = 'pypath',
        level = loglevel,
        wrap = False,
    )

# include info at the beginning of the log:
info(0)
