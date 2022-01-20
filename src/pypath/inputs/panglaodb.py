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
#           Olga Ivanova
#           Sebastian Lobentanzer
#           Ahmet Rifaioglu
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

import csv
import collections

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.settings as settings
import pypath.share.common as common
import pypath.utils.mapping as mapping


def panglaodb_raw():
    """
    Cell type marker genes from PanglaoDB (https://panglaodb.se/index.html).
    """

    url = urls.urls['panglaodb']['url']
    headers = [settings.get('user_agent')]
    c = curl.Curl(url, silent = False, large = True, req_headers = headers)

    result = list(csv.DictReader(c.result, delimiter = '\t'))

    return result


def panglaodb_annotations():
    """
    Cell type marker genes from PanglaoDB (https://panglaodb.se/index.html).
    """

    to_float = common.float_or_nan

    record = collections.namedtuple(
        'PanglaodbAnnotation',
        (
            'canonical_marker',
            'cell_type',
            'organ',
            'germ_layer',
            'entity_type',
            'human',
            'mouse',
            'human_sensitivity',
            'mouse_sensitivity',
            'human_specificity',
            'mouse_specificity',
            'ubiquitiousness',
        )
    )

    result = collections.defaultdict(set)

    raw = panglaodb_raw()

    for r in raw:

        uniprots = mapping.map_name(
            r['official gene symbol'],
            'genesymbol',
            'uniprot',
        )

        for u in uniprots:

            result[u].add(
                record(
                    canonical_marker = r['canonical marker'] == '1',
                    cell_type = r['cell type'],
                    organ = r['organ'],
                    germ_layer = r['germ layer'],
                    entity_type = 'protein', # TODO: handle other entity types
                    human = 'Hs' in r['species'],
                    mouse = 'Mm' in r['species'],
                    human_sensitivity = to_float(r['sensitivity_human']),
                    mouse_sensitivity = to_float(r['sensitivity_mouse']),
                    human_specificity = to_float(r['specificity_human']),
                    mouse_specificity = to_float(r['specificity_mouse']),
                    ubiquitiousness = to_float(r['ubiquitousness index']),
                )
            )

    return result
