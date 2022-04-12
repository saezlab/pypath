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

import collections

import bs4

import pypath.resources.urls as urls
import pypath.share.curl as curl


"""
Note: find other Ensembl related functions in ``inputs.biomart``.
"""


def ensembl_organisms():
    """
    List of organisms in Ensembl with various taxon IDs and metadata about
    related Ensembl database contents.

    Returns:
        List of named tuples.
    """

    record = None
    result = []
    url = urls.urls['ensembl']['species']
    c = curl.Curl(url)
    soup = bs4.BeautifulSoup(c.result, 'html.parser')

    for r in soup.find('table').find_all('tr'):

        if not record:

            record = collections.namedtuple(
                'EnsemblOrganism',
                [c.text.lower().replace(' ', '_') for c in r] +
                ['ensembl_name']
            )

            continue

        r = list(r)

        result.append(
            record(
                *(
                    int(c.text) if i == 2 else c.text
                    for i, c in enumerate(r)
                ),
                (
                    # Mus musculus -> mmusculus
                    lambda x: ''.join([xx[0] for xx in x[:-1]] + [x[-1]])
                )(
                    r[1].text.lower().split()
                )
            )
        )

    return result
