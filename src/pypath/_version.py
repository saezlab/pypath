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

__all__ = [
    '_get_version',
    '__version__',
    '__release__',
    '__author__',
    '_MINOR',
    '_MAJOR',
    '_MICRO',
]

import os

def _get_version():
    """
    Returns the version number.

    Reads file named ``__version__`` and returns it as a tuple of
    integers.

    :return:
        (*tuple*) -- Three element tuple containing the major, minor and
        micro version numbers.

    **Example:**
        >>> # If ``__version__`` contains ``1.11.0``:
        >>> _get_version()
        (1, 11, 0)
    """

    dir_root = os.path.abspath(os.path.dirname(__file__))

    with open(os.path.join(dir_root, '__version__'), 'r') as v:
        
        return tuple([int(i) for i in v.read().strip().split('.')])


_MAJOR, _MINOR, _MICRO = _get_version()
__version__ = '%d.%d.%d' % (_MAJOR, _MINOR, _MICRO)
__release__ = '%d.%d' % (_MAJOR, _MINOR) # XXX: Not used
__author__ = (
    u'Dénes Türei, '
    u'Nicolàs Palacio, '
    u'Olga Ivanova. \n'
    u'European Molecular Biology Laboratory, Heidelberg Germany\n'
    u'European Bioinformatics Institute, Hinxton UK\n'
    u'University Hospital RWTH, Aachen Germany\n'
    u'University Hospital Heidelberg Germany'
)
