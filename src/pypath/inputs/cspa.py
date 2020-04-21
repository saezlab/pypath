#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2020
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#                  Nicolàs Palacio
#                  Olga Ivanova
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

import collections

import pypath.share.curl as curl
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
        ]
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

    return result


def cspa_cell_type_annotations():
