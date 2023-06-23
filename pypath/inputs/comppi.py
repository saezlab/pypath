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

import itertools
import collections

import pypath.share.curl as curl
import pypath.share.session as session
import pypath.share.settings as settings
import pypath.resources.urls as urls
import pypath.utils.mapping as mapping
import pypath.utils.taxonomy as taxonomy

_logger = session.Logger(name = 'comppi_input')
_log = _logger._log


def comppi_interaction_locations(organism = 9606):
    """
    Downloads and preprocesses protein interaction and cellular compartment
    association data from the ComPPI database.
    This data provides scores for occurrence of protein-protein interactions
    in various compartments.
    """

    ComppiLocation = collections.namedtuple(
        'ComppiLocation',
        [
            'location',
            'score',
        ],
    )

    ComppiInteraction = collections.namedtuple(
        'ComppiInteraction',
        [
            'id_a',
            'id_b',
            'loc_a',
            'loc_b',
        ],
    )

    def process_locations(loc):

        return tuple(
            ComppiLocation(location = llloc[0], score = float(llloc[1]))
            for llloc in
            (lloc.split(':') for lloc in loc.split('|'))
        )


    organisms = {
        9606: 0,
        7227: 1,
        6239: 2,
        4932: 3,
    }

    ncbi_tax_id = taxonomy.ensure_ncbi_tax_id(organism)

    if ncbi_tax_id not in organisms:

        raise ValueError(
            'Can not recognize organism: `%s`. '
            'Available organisms are human (9606), Drosophila (7227), '
            'C. elegans (6239) and S. cerevisiae (4932).' % str(organism)
        )

    url = urls.urls['comppi']['url']
    headers = [settings.get('user_agent')]

    # obtaining cookie
    c = curl.Curl(url, cache = False, req_headers = headers)
    cookie = ';'.join([
        h.decode().split(':')[1].split(';')[0].strip()
        for h in c.resp_headers
        if h.startswith(b'Set-Cookie')
    ])
    cookie_hdr = ['Cookie: %s' % cookie]
    _log('Cookie from ComPPI: %s' % cookie)

    # performing the download
    post = {
        'fDlSet': 'comp',
        'fDlSpec': '%u' % ncbi_tax_id,
        'fDlMLoc': 'all',
        'fDlSubmit': 'Download'
    }
    c = curl.Curl(
        url,
        req_headers = headers + cookie_hdr,
        post = post,
        large = True,
        silent = False,
        compr = 'gz',
    )

    _ = next(c.result)

    for l in c.result:

        l = l.strip('\r\n').split('\t')

        organism_a = int(l[7])
        organism_b = int(l[15])

        if organism and (organism_a != organism or organism_b != organism):
            continue

        for uniprot1, uniprot2 in itertools.product(
            mapping.map_name(l[0], 'uniprot', 'uniprot'),
            mapping.map_name(l[8], 'uniprot', 'uniprot'),
        ):

            yield ComppiInteraction(
                id_a = uniprot1,
                id_b = uniprot2,
                loc_a = process_locations(l[2]),
                loc_b = process_locations(l[10]),
            )


def comppi_locations(organism = 9606, score_threshold = .7):

    result = collections.defaultdict(set)

    for iloc in comppi_interaction_locations(organism = organism):

        for label in ('a', 'b'):

            for loc in getattr(iloc, 'loc_%s' % label):

                if loc.location == 'N/A' or loc.score < score_threshold:

                    continue

                result[getattr(iloc, 'id_%s' % label)].add(loc)

    return dict(result)
