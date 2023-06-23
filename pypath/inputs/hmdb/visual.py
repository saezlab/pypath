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
Visualisation services from the Human Metabolome Database (HMDB).
"""

import webbrowser

import pypath.resources.urls as urls
import pypath.share.curl as curl


def _svg_url(hmdb_id: str) -> str:
    """
    Returns the SVG URL for the given HMDB ID.
    """

    return urls.urls['hmdb']['svg'] % hmdb_id


def show_structure(hmdb_id: str) -> None:
    """
    Show the 2D structure of the given HMDB ID in a web browser.

    Args:
        hmdb_id:
            Identifier of a metabolite.
    """

    url = _svg_url(hmdb_id)
    webbrowser.open(url)


def structure_svg(hmdb_id: str, path: str | None = None) -> str:
    """
    Retrieve the picture of the 2D structure of a metabolite in SVG format.

    Args:
        hmdb_id:
            Identifier of a metabolite.
        path:
            Save the SVG file to this path (optional).

    Returns:
        The SVG file as a string.
    """

    url = _svg_url(hmdb_id)

    c = curl.Curl(url, large = False, silent = True)

    if path:

        with open(path, 'w') as fp:

            fp.write(c.result)

    return c.result
