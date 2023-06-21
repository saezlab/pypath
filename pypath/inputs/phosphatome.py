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
import pypath.utils.mapping as mapping
import pypath.inputs.common as inputs_common
import pypath.inputs.science as science_input


def phosphatome_annotations():
    """
    Downloads the list of phosphatases from Chen et al, Science Signaling
    (2017) Table S1.
    """

    PhosphatomeAnnotation = collections.namedtuple(
        'PhosphatomeAnnotation',
        [
            'fold',
            'family',
            'subfamily',
            'has_protein_substrates',
            'has_non_protein_substrates',
            'has_catalytic_activity',
        ],
    )

    url = urls.urls['phosphatome']['url']
    path = science_input.science_download(url = url)
    c = curl.FileOpener(
        path,
        compr = 'zip',
        files_needed = ['aag1796_Tables S1 to S23.xlsx'],
        large = True,
        default_mode = 'rb',
    )
    tbl = inputs_common.read_xls(c.result['aag1796_Tables S1 to S23.xlsx'])
    result = []

    result = collections.defaultdict(set)

    for rec in tbl[2:]:

        uniprots = mapping.map_name(rec[0], 'genesymbol', 'uniprot')

        for uniprot in uniprots:

            result[uniprot].add(
                PhosphatomeAnnotation(
                    fold = rec[2],
                    family = rec[3],
                    subfamily = rec[4],
                    has_protein_substrates = rec[21].strip().lower() == 'yes',
                    has_non_protein_substrates = (
                        rec[22].strip().lower() == 'yes'
                    ),
                    has_catalytic_activity = rec[23].strip().lower() == 'yes',
                )
            )

    return dict(result)
