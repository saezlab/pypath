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

import re
import collections

import bs4

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.utils.taxonomy as taxonomy


def pdzbase_interactions():
    """
    Downloads data from PDZbase. Parses data from the HTML tables.

    Returns
        List of named tuples with interaction data.
    """

    PDZbaseInteraction = collections.namedtuple(
        'PDZbaseInteraction',
        [
            'uniprot_pdz',
            'isoform_pdz',
            'uniprot_ligand',
            'isoform_ligand',
            'genesymbol_pdz',
            'genesymbol_ligand',
            'pdz_domain',
            'organism',
            'pubmed',
        ],
    )

    # UniProt ID with isoform e.g. O14754-1
    reupi = re.compile(r'([\w]{6,10})(?:-([0-9]{1,2}))?')

    url = urls.urls['pdzbase']['url_rescued']
    c = curl.Curl(url, silent = False)
    data = c.result
    soup = bs4.BeautifulSoup(data, 'html.parser')
    rows = (
        soup.find_all('table')[3].find('table').find('table').find_all('tr')
    )
    result = []

    del rows[0]

    for r in rows:

        r = [c.text.strip() for c in r.find_all('td')]

        uniprot_pdz, isoform_pdz = reupi.match(r[1]).groups()
        uniprot_ligand, isoform_ligand = reupi.match(r[4]).groups()

        result.append(
            PDZbaseInteraction(
                uniprot_pdz = uniprot_pdz,
                isoform_pdz = int(isoform_pdz) if isoform_pdz else 1,
                uniprot_ligand = uniprot_ligand,
                isoform_ligand = int(isoform_ligand) if isoform_ligand else 1,
                genesymbol_pdz = r[0],
                genesymbol_ligand = r[3],
                pdz_domain = int(r[2]),
                organism = taxonomy.ensure_ncbi_tax_id(r[5]),
                pubmed = int(r[6]),
            )
        )

    return result

