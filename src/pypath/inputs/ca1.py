#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2021
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

import collections

import pypath.resources.urls as urls
import pypath.share.curl as curl


def ca1_interactions():
    """
    Downloads and processes the CA1 signaling network (Ma\'ayan 2005).
    Returns list of interactions.
    """

    Ca1Interaction = collections.namedtuple(
        'Ca1Interaction',
        (
            'source_label',
            'source_uniprot',
            'source_uniprot_mouse',
            'source_function',
            'source_location',
            'target_label',
            'target_uniprot',
            'target_uniprot_mouse',
            'target_function',
            'target_location',
            'effect',
            'interaction_type',
            'pmid',
        ),
    )

    url = urls.urls['ca1']['url']
    c = curl.Curl(url, silent = False, files_needed = ['S1.txt'])
    data = c.result
    result = []

    for l in data['S1.txt'].split('\n')[1:]:

        l = l.strip().split()

        if len(l) == 13:

            result.append(Ca1Interaction(*l))

    return result
