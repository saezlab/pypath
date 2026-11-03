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

import csv
import collections

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.utils.mapping as mapping
import pypath.share.common as common


LrdbAnnotation = collections.namedtuple(
    'LrdbAnnotation',
    [
        'role',
        'cell_type',
        'sources',
        'references',
    ],
)


def lrdb_interactions():

    resource_names = {
        'reactome': 'Reactome',
        'fantom5': 'Fantom5',
        'IUPHAR': 'Guide2Pharma',
        'uniprot': 'UniProt',
    }

    def remove(lst, to_remove):
        to_remove = common.to_set(to_remove)

        return [
            it
            for it in lst
            if it not in to_remove
        ]

    LrdbRecord = collections.namedtuple(
        'LrdbRecord',
        [
            'ligand_genesymbol',
            'receptor_genesymbol',
            'sources',
            'references',
            'ligand_cells',
            'receptor_cells',
        ]
    )

    url = urls.urls['lrdb']['url']

    c = curl.Curl(url, silent = False, large = True)

    reader = csv.DictReader(c.result, delimiter = '\t')

    result = []

    for rec in reader:

        result.append(
            LrdbRecord(
                ligand_genesymbol = rec['ligand'],
                receptor_genesymbol = rec['receptor'],
                sources = [
                    resource_names[src] if src in resource_names else src
                    for src in
                    remove(
                        rec['source'].split(','),
                        {'literature', ''},
                    )
                ],
                references = remove(rec['PMIDs'].split(','), ''),
                ligand_cells = remove(rec['cells.L'].split(','), ''),
                receptor_cells = remove(rec['cells.R'].split(','), ''),
            )
        )

    return result


def lrdb_annotations():

    result = collections.defaultdict(set)

    lrdb = lrdb_interactions()

    for rec in lrdb:

        for role in ('ligand', 'receptor'):

            uniprots = mapping.map_name(
                getattr(rec, '%s_genesymbol' % role),
                'genesymbol',
                'uniprot',
            )

            for uniprot in uniprots:

                cell_types = getattr(rec, '%s_cells' % role) or (None,)

                for cell_type in cell_types:

                    cell_type = (
                        'T lymphocyte'
                            if cell_type == 'tymphocyte' else
                        cell_type.replace('cells', 'cell')
                            if cell_type else
                        None
                    )

                    result[uniprot].add(
                        LrdbAnnotation(
                            role = role,
                            cell_type = cell_type,
                            sources = tuple(sorted(rec.sources)),
                            references = tuple(sorted(rec.references)),
                        )
                    )

    return dict(result)
