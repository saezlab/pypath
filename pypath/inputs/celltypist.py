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

import pypath.inputs.common as inputs_common
import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.utils.mapping as mapping


def celltypist_annotations():
    """
    Immune cell type marker genes from Celltypist
    (https://github.com/Teichlab/celltypist_wiki).
    """

    recomma = re.compile(r'\s?,\s?')

    def multi_value_field(value):

        return tuple(sorted(recomma.split(value.strip())))

    record = collections.namedtuple(
        'CelltypistAnnotation',
        (
            'cell_type',
            'cell_subtype',
            'cell_ontology',
            'marker_type',
            'tissues',
            'datasets',
        )
    )

    result = collections.defaultdict(set)

    url = urls.urls['celltypist']['url']

    c = curl.Curl(url, silent = False, large = True)

    xls = c.fileobj
    xlsfile = xls.name
    xls.close()
    tbl = inputs_common.read_xls(xlsfile)[1:]

    marker_columns = ((6, 'curated_marker'), (7, 'celltypist_model'))

    for r in tbl:

        for col, marker_type in marker_columns:

            genesymbols = recomma.split(r[col].strip())
            uniprots = mapping.map_names(genesymbols, 'genesymbol', 'uniprot')

            annot = record(
                cell_type = r[0],
                cell_subtype = r[1],
                cell_ontology = r[3],
                marker_type = marker_type,
                tissues = multi_value_field(r[4]),
                datasets = multi_value_field(r[5]),
            )

            for u in uniprots:

                result[u].add(annot)

    return dict(result)
