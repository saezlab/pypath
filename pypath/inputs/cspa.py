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

from future.utils import iteritems

import collections

import pypath.share.curl as curl
import pypath.share.common as common
import pypath.resources.urls as urls
import pypath.inputs.common as inputs_common
import pypath.utils.mapping as mapping
import pypath.utils.taxonomy as taxonomy


def cspa_annotations(organism = 9606):


    CspaAnnotation = collections.namedtuple(
        'CspaAnnotation',
        [
            'high_confidence',
            'n_cell_types',
            'tm',
            'gpi',
            'uniprot_cell_surface',
        ],
    )


    sheets = {
        'Human': 'Table A',
        'Mouse': 'Table B',
    }

    str_organism = taxonomy.taxids[organism].capitalize()

    url = urls.urls['cspa']['url_s2']
    c = curl.Curl(url, large = True, silent = False)
    xlsname = c.fname
    del(c)
    raw = inputs_common.read_xls(xlsname, sheets[str_organism])[1:]

    result = collections.defaultdict(set)

    for row in raw:

        for uniprot in mapping.map_name(row[1], 'uniprot', 'uniprot'):

            result[uniprot].add(
                CspaAnnotation(
                    high_confidence = 'high confidence' in row[2],
                    n_cell_types = int(float(row[9])),
                    tm = int(float(row[11])),
                    gpi = int(float(row[12])),
                    uniprot_cell_surface = row[13] == 'yes',
                )
            )

    return dict(result)


def cspa_cell_types(organism = 9606):

    sheets = {
        'Human': 'Table_E',
        'Mouse': 'Table_F',
    }

    str_organism = taxonomy.taxids[organism].capitalize()

    url = urls.urls['cspa']['url_s1']
    c = curl.Curl(url, large = True, silent = False)
    xlsname = c.fname
    del(c)
    raw = inputs_common.read_xls(xlsname, sheets[str_organism])

    result = collections.defaultdict(lambda: collections.defaultdict(dict))

    cell_types = raw[0][1:]

    for row in raw[1:]:

        for uniprot in mapping.map_name(row[0], 'uniprot', 'uniprot'):

            for col, cell_type in enumerate(cell_types):

                value = row[col + 1]

                result[cell_type][uniprot] = (
                    float(value)
                        if common.is_float(value) else
                    None
                )

    return dict((k, dict(v)) for k, v in iteritems(result))


def cspa_cell_type_annotations(organism = 9606):


    CspaCellType = collections.namedtuple(
        'CspaCellType',
        [
            'cell_type',
            'value',
        ],
    )


    cell_type_data = cspa_cell_types(organism = organism)


    result = collections.defaultdict(set)

    for cell_type, data in iteritems(cell_type_data):

        for uniprot, value in iteritems(data):

            if value:

                result[uniprot].add(
                    CspaCellType(
                        cell_type = cell_type,
                        value = value,
                    )
                )

    return dict(result)
