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

import pypath.resources.urls as urls
import pypath.share.curl as curl


def lincs_compounds():
    """
    The returned dict has names, brand names or company specific IDs of
    compounds as keys, and tuples of PubChem, ChEMBL, ChEBI, InChi,
    InChi Key, SMILES and LINCS as values.
    """

    LincsCompound = collections.namedtuple(
        'LincsCompound',
        (
            'lincs',
            'chembl',
            'chebi',
            'inchi',
            'inchi_key',
            'smiles',
            'alternatives',
        )
    )

    c = curl.Curl(urls.urls['lincs-compounds']['url'], silent = False)

    return dict(
        [
            (key, pair[1])
            for pair in [
                (
                    [
                        it for sl in
                        [
                            filter(
                                lambda z: len(z) > 0,
                                y.split(';')
                            )
                            for y in x[1:4]
                            if len(y) > 0
                        ]
                        for it in sl
                    ],
                    LincsCompound(
                        lincs = x[4],
                        chembl = 'CHEMBL%s' % x[7] if x[7] else None,
                        chebi = 'CHEBI:%s' % x[8] if x[8] else None,
                        inchi = x[9],
                        inchi_key = x[10],
                        smiles = x[11],
                        alternatives = x[3],
                    )
                )
                for x in [
                    [b.strip() for b in a.split('\t')]
                    for a in ''.join([
                        s.replace(',', '\t')
                            if i % 2 == 0 else
                        s.replace('\n', '')
                        for i, s in enumerate(c.result.split('"'))
                    ]).split('\n')[1:]
                    if len(a) > 0
                ]
            ]
            for key in pair[0]
        ]
    )