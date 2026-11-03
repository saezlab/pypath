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

from numbers import Number
from typing import List

import collections

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.progress as progress


def intact_interactions(
        miscore: Number = .6,
        organism: int = 9606,
        complex_expansion: bool = False,
        only_proteins: bool = False,
        only_ids: bool = False,
    ) -> List[tuple]:
    """
    only_proteins : bool
        Keep only records of protein-protein interactions.
    only_ids : bool
        Load only the identifiers of interacting pairs
        (smaller memory footprint).
    """

    id_types = {
        'uniprotkb': 'uniprot',
    }

    IntactInteraction = collections.namedtuple(
        'IntactInteraction',
        (
            'id_a',
            'id_b',
            'id_type_a',
            'id_type_b',
            'pubmeds',
            'methods',
            'interaction_types',
            'mi_score',
            'isoform_a',
            'isoform_b',
        ),
    )
    IntactInteraction.__new__.__defaults__ = (None,) * 7


    def get_id_type(field):

        id_type = None if field == '-' else field.split(':')[0]

        return id_types[id_type] if id_type in id_types else id_type


    def get_uniprot_id(field):

        uniprot, isoform = _try_isoform(
            field.split(':')[1].replace('"', '')
        )

        uniprot = uniprot.split('-')[0]

        return uniprot, isoform

    def get_ebi_id(field):

        if field == '-':

            return None, None
        
        else:

            partner_id=field.split(':')[1]

            return partner_id, None


    def get_taxon(field):

        return (
            0
                if field == '-' else
            field.split('|')[0].split(':')[1].split('(')[0]
        )


    results = []
    url = urls.urls['intact']['mitab']

    if type(organism) is int:
        organism = '%u' % organism

    c = curl.Curl(
        url,
        silent = False,
        large = True,
        files_needed = ['intact.txt'],
        slow = True,
    )

    data = c.result['intact.txt']
    size = c.sizes['intact.txt']
    prg = progress.Progress(size, 'Reading IntAct MI-tab file', 99)

    for lnum, l in enumerate(data):

        prg.step(len(l))

        if lnum == 0:

            continue

        l = l.strip('\n\r ').split('\t')

        taxon_a = get_taxon(l[9])
        taxon_b = get_taxon(l[10])

        if (
            (
                organism is None or (
                    taxon_a == organism and
                    taxon_b == organism
                )
            ) and (
                complex_expansion or
                'expansion' not in l[15]
            )
        ):

            # finding mi-score and author
            sc = 0
            au = '0'

            for s in l[14].split('|'):

                if s.startswith('intact-miscore'):

                    sc = float(s.split(':')[1])

                if s.startswith('author'):

                    au = len(s.split(':')[1])

            # filtering for mi-score
            if sc < miscore:

                continue

            id_type_a = get_id_type(l[0])
            id_type_b = get_id_type(l[1])

            if (
                only_proteins and not (
                    id_type_a == 'uniprot' and
                    id_type_b == 'uniprot'
                )
            ):

                continue

            id_a, isoform_a = (
                get_uniprot_id(l[0])
                    if id_type_a == 'uniprot' else
                get_ebi_id(l[0])
            )

            id_b, isoform_b = (
                get_uniprot_id(l[1])
                    if id_type_b == 'uniprot' else
                get_ebi_id(l[1])
            )

            # key = tuple(sorted((id_a, id_b)))

            pubmeds = set(
                ref[1] for ref in (
                    ref.split(':')
                    for ref in l[8].split('|')
                )
                if ref[0] == 'pubmed'
            )
            methods = set(
                met.split('(')[1].strip(')"')
                for met in  l[6].split('|')
            )

            interaction_types= set(
                int_type.split('(')[1].strip(')"')
                for int_type in  l[11].split('|')
            )

            results.append(
                IntactInteraction(
                    id_a = id_a,
                    id_b = id_b,
                    id_type_a = id_type_a,
                    id_type_b = id_type_b,
                    pubmeds = pubmeds,
                    methods = methods,
                    interaction_types = interaction_types,
                    mi_score = sc,
                    isoform_a = isoform_a,
                    isoform_b = isoform_b,
                )
            )

    prg.terminate()

    return results


def _try_isoform(name):

    name = name.split('-')

    if len(name) > 1 and name[1].isdigit():

        isoform = int(name[1])
        main = name[0]

    else:

        main = '-'.join(name)
        isoform = None

    return main, isoform
