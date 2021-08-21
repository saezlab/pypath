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


def acsn_interactions(keep_in_complex_interactions = True):
    """
    Processes ACSN data from local file.
    Returns list of interactions.

    @keep_in_complex_interactions : bool
        Whether to include interactions from complex expansion.
    """

    AcsnInteraction = collections.namedtuple(
        'AcsnInteraction',
        (
            'partner_a',
            'partner_b',
            'mechanism',
            'references',
        ),
    )

    names_url = urls.urls['acsn']['names']
    ppi_url = urls.urls['acsn']['ppi']
    names_c = curl.Curl(names_url, silent = False, large = True)
    ppi_c = curl.Curl(ppi_url, silent = False, large = True)

    names = {}
    interactions = []

    for l in names_c.result:

        l = l.strip().split('\t')
        names[l[0]] = l[2:]

    _ = next(ppi_c.result)

    for l in ppi_c.result:

        l = l.strip().split('\t')

        if l[0] in names:

            for a in names[l[0]]:

                if l[2] in names:

                    for b in names[l[2]]:

                        if keep_in_complex_interactions:

                            if 'PROTEIN_INTERACTION' in l[1]:

                                l[1].replace(
                                    'COMPLEX_EXPANSION',
                                    'IN_COMPLEX_INTERACTION'
                                )

                        interactions.append(
                            AcsnInteraction(
                                partner_a = a,
                                partner_b = b,
                                mechanism = l[1],
                                references = l[3]
                            )
                        )

    return interactions
