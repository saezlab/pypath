#!/usr/bin/env python

#
# This file is part of the `omnipath_utils` Python module
#
# Copyright 2026
# Heidelberg University Hospital
#
# File author(s): Denes Turei (turei.denes@gmail.com)
#
# Distributed under the GPL-3.0-or-later license
# See the file `LICENSE` or read a copy at
# https://www.gnu.org/licenses/gpl-3.0.txt
#

"""Package metadata (version, authors, etc)."""

__all__ = ['__version__', '__author__', '__license__']

import importlib.metadata

_FALLBACK_VERSION = '0.1.1'

try:
    __version__ = importlib.metadata.version('omnipath_utils')
except importlib.metadata.PackageNotFoundError:
    # Package not installed (e.g. running from source checkout)
    __version__ = _FALLBACK_VERSION

__author__ = 'Denes Turei'
__license__ = 'GPL-3.0-or-later'
