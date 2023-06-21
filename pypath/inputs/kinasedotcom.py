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

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.utils.mapping as mapping
import pypath.inputs.common as inputs_common


def kinasedotcom_annotations():
    """
    Downloads and processes kinase annotations from kinase.com.
    """

    KinasedotcomAnnotation = collections.namedtuple(
        'KinasedotcomAnnotation',
        ['group', 'family', 'subfamily']
    )
    KinasedotcomAnnotation.__new__.__defaults__ = (None,)


    def add_record(uniprot, rec, offset = 2):

        if rec[offset].strip():
            result[uniprot].add(
                KinasedotcomAnnotation(
                    group = rec[offset].strip(),
                    family = rec[offset + 1].strip(),
                    subfamily = rec[offset + 2].strip() or None,
                )
            )


    url = urls.urls['kinome']['rescued']
    c = curl.Curl(url, large = True, silent = False)
    xlsf = c.fileobj
    xlsname = xlsf.name
    xlsf.close()
    tbl = inputs_common.read_xls(xlsname)

    result = collections.defaultdict(set)

    for rec in tbl:

        uniprots = mapping.map_name(rec[23].strip(), 'genesymbol', 'uniprot')

        for uniprot in uniprots:

            add_record(uniprot, rec)

            if rec[12].strip():

                add_record(uniprot, rec, offset = 12)

    return dict(result)
