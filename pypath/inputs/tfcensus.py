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

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.utils.mapping as mapping


def tfcensus_annotations(only_classes = None):
    """
    Downloads and processes the list of all known transcription factors from
    TF census (Vaquerizas 2009). This resource is human only.
    Returns dict with UniProt IDs as keys and TF annotations as values.
    """

    TfcensusAnnotation = collections.namedtuple(
        'TfcensusAnnotation',
        ['tfcensus_class', 'tissue'],
    )


    result = collections.defaultdict(set)

    reensg = re.compile(r'ENSG[0-9]{11}')
    url = urls.urls['vaquerizas2009']['url']
    c = curl.Curl(url, silent = False, large = True)

    header = True

    for l in c.result:

        l = l.split('\t')

        if header:

            if l[0] == 'Class':

                header = False

            continue

        ensg = reensg.findall(l[1])
        hgnc = l[5]
        tfcensus_class = l[0]
        tissues = l[6].strip() if len(l) > 6 else ''
        tissues = tissues.split(';') if tissues else [None]

        uniprots = mapping.map_names(ensg, 'ensembl', 'uniprot')
        uniprots.update(mapping.map_name(hgnc, 'genesymbol', 'uniprot'))

        for uniprot in uniprots:

            for tissue in tissues:

                result[uniprot].add(
                    TfcensusAnnotation(
                        tfcensus_class = tfcensus_class,
                        tissue = tissue,
                    )
                )

    return dict(result)