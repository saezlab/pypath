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

from future.utils import iteritems

import collections

import pypath.share.curl as curl
import pypath.share.settings as settings
import pypath.resources.urls as urls
import pypath.utils.taxonomy as taxonomy


def get_exocarta(organism = 9606, types = None):
    """
    :param set types:
        Molecule types to retrieve. Possible values: `protein`, `mrna`.
    """

    return _get_exocarta_vesiclepedia(
        database = 'exocarta',
        organism = organism,
        types = types,
    )

def get_vesiclepedia(organism = 9606, types = None):
    """
    :param set types:
        Molecule types to retrieve. Possible values: `protein`, `mrna`.
    """

    return _get_exocarta_vesiclepedia(
        database = 'vesiclepedia',
        organism = organism,
        types = types,
    )

def _get_exocarta_vesiclepedia(
        database = 'exocarta',
        organism = 9606,
        types = None
    ):
    """
    :param str database:
        Which database to download: ExoCarta or Vesiclepedia.
    :param set types:
        Molecule types to retrieve. Possible values: `protein`, `mrna`.
    """

    database = database.lower()

    types = types or {'protein'}

    organism = taxonomy.phosphoelm_taxids[organism]

    taxid_rev = dict((v, k) for k, v in iteritems(taxonomy.phosphoelm_taxids))

    # collecting the references
    url_s = urls.urls[database]['url_study']
    c = curl.Curl(url_s, large = True, silent = False)
    _ = next(c.result)

    studies = {}

    for s in c.result:
        s = s.split('\t')

        organisms = tuple(
            taxid_rev[t.strip()]
            for t in s[2].split('|')
            if t.strip() in taxid_rev
        )

        if not organisms:
            continue

        stud = (
            s[1] if s[1] != '0' else None, # PubMed ID
            organisms, # organism
            s[4], # sample source (cell type, tissue)
        )

        if database == 'vesiclepedia':
            vtype = s[11].strip()

            stud += (
                tuple(vtype.split('/')) if vtype else (),
            )

        studies[int(s[0])] = tuple(stud)

    # processing proteins
    url_p = urls.urls[database]['url_protein']
    c = curl.Curl(url_p, large = True, silent = False, slow = True)
    _ = next(c.result)

    for s in c.result:
        s = s.split('\t')

        if s[4] != organism or s[1] not in types:
            continue

        yield (
            s[2], # Entrez ID
            s[3], # Gene Symbol
            taxid_rev[s[4]], # NCBI Taxonomy ID
            studies[int(s[5])], # study reference
        )
