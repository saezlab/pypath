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

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.inputs.common as inputs_common


def mirecords_interactions():
    """
    Retrieves literature curated miRNA-target gene interactions from
    miRecords (c1.accurascience.com/miRecords/).
    """

    MirecordsInteraction = collections.namedtuple(
        'MirecordsInteraction',
        (
            'mirna_name',
            'target_refseq',
            'target_genesymbol',
            'mirna_organism',
            'target_organism',
            'pmid',
        ),
    )

    url = urls.urls['mirecords']['url']
    c = curl.Curl(url, silent = False, large = True)

    tbl = inputs_common.read_xls(c.fileobj.name)

    c.close()

    return [
        MirecordsInteraction(
            mirna_name = l[6],
            target_refseq = l[3],
            target_genesymbol = l[2],
            mirna_organism = l[1],
            target_organism = l[5],
            pmid = l[0].split('.')[0],
        )
        for l in
        ([f.strip() for f in ll] for ll in tbl[1:])
    ]
