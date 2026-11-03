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

import collections
import itertools

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.common as common
import pypath.utils.mapping as mapping
import pypath.utils.taxonomy as taxonomy


def hippie_interactions(
        score_threshold = .75,
        only_human = False,
        only_sources = None,
        only_methods = None,
        methods = False,
        sources = False,
        references = True,
        organisms = False,
    ):

    only_sources = common.to_set(only_sources)
    only_methods = common.to_set(only_methods)

    HippieInteraction = collections.namedtuple(
        'HippieInteraction',
        [
            'id_a',
            'id_b',
            'score',
            'methods',
            'references',
            'sources',
            'organisms',
        ],
    )

    tps = lambda i: tuple(sorted(i))

    url = urls.urls['hippie']['url']
    c = curl.Curl(url, large = True, silent = False)

    result = set()

    for i, l in enumerate(c.result):

        l = l.strip('\r\n').split('\t')

        score = float(l[4])

        if score < score_threshold:

            continue

        ids_a_1 = mapping.map_name(l[0], 'uniprot-entry', 'uniprot')
        ids_a_2 = mapping.map_name(l[1], 'entrez', 'uniprot')
        ids_b_1 = mapping.map_name(l[2], 'uniprot-entry', 'uniprot')
        ids_b_2 = mapping.map_name(l[3], 'entrez', 'uniprot')

        for id_a, id_b in itertools.product(
            ids_a_1 | ids_a_2,
            ids_b_1 | ids_b_2
        ):

            details = dict(
                (
                    dd[0],
                    set(dd[1].split(',')),
                )
                for dd in
                (d.split(':') for d in l[5].split(';'))
            )

            _sources = details['sources'] if 'sources' in details else set()
            experiments = (
                details['experiments'] if 'experiments' in details else set()
            )

            if not all((
                not only_methods or experiments & only_methods,
                not only_methods or _sources & only_sources,
            )):

                continue

            _organisms = {9606}

            if 'species' in details:

                names = {
                    spec.split('(')[0].strip()
                    for spec in details['species']
                }
                _organisms = {
                    taxonomy.ensure_ncbi_tax_id(name)
                    for name in names
                }
                _organisms.discard(None)

                if only_human and 9606 not in _organisms:

                    continue

            result.add(
                HippieInteraction(
                    id_a = id_a,
                    id_b = id_b,
                    score = score,
                    methods = tps(experiments) if methods else None,
                    references = (
                        tps(details['pmids']) if references else None
                    ),
                    sources = tps(_sources) if sources else None,
                    organisms = tps(_organisms) if organisms else None,
                )
            )

    return list(result)
