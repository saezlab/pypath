#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2021
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: Dénes Türei (turei.denes@gmail.com)
#           Nicolàs Palacio
#           Olga Ivanova
#           Sebastian Lobentanzer
#           Ahmet Rifaioglu
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

import pypath.resources.urls as urls
import pypath.share.curl as curl


def get_hsn():
    """
    Downloads and processes HumanSignalingNetwork version 6
    (published 2014 Jan by Edwin Wang).
    Returns list of interactions.
    """

    url = urls.urls['hsn']['url']
    c = curl.Curl(url, silent = False, large = True)
    data = c.result
    data = [r.strip().split(',') for r in data if r][1:]

    return data


def wang_interactions():
    """
    Downloads and processes Wang Lab HumanSignalingNetwork.
    Returns list of interactions as tuples of source, target and effect.
    """

    url = urls.urls['wang']['rescued']
    c = curl.Curl(url, silent = False)
    data = c.result
    data = data.split('\n')
    effects = []
    nodes = {}
    reading_nodes = False
    reading_edges = False

    for l in data:

        if len(l.strip()) == 0:
            reading_nodes = False
            reading_edges = False

        l = l.split(',')

        if reading_nodes:
            nodes[l[0]] = l[1]

        if reading_edges:
            effects.append([nodes[l[0]], nodes[l[1]], l[2]])

        if l[0].startswith('Node'):
            reading_nodes = True

        if l[0].startswith('From'):
            reading_nodes = False
            reading_edges = True

    return effects
