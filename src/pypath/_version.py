#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright (c) 2014-2017 - EMBL-EBI
#
#  File author(s): Dénes Türei (denes@ebi.ac.uk)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://www.ebi.ac.uk/~denes
#

import os

def _get_version():
    """Returns the version number.

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

    # XXX: Why all caps? this is not a global variable
    ROOT = os.path.abspath(os.path.dirname(__file__))

    with open(os.path.join(ROOT, '__version__'), 'r') as v:
        return tuple([int(i) for i in v.read().strip().split('.')])

_MAJOR, _MINOR, _MICRO = _get_version()
__version__ = '%d.%d.%d' % (_MAJOR, _MINOR, _MICRO)
__release__ = '%d.%d' % (_MAJOR, _MINOR) # XXX: Not used
__author__ = u'Dénes Türei'
