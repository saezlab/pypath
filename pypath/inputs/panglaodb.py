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
            'ncbi_tax_id',
        )
    )

    result = collections.defaultdict(set)

    raw = panglaodb_raw()

    for r in raw:

        human = 'Hs' in r['species'],
        mouse = 'Mm' in r['species'],

        entity_type = 'protein'
        uniprots = mapping.map_name(
            r['official gene symbol'],
            'genesymbol',
            'uniprot',
        )
        this_ncbi_tax_id = 9606

        if not uniprots and human:

            uniprots = mapping.map_name(
                r['official gene symbol'],
                'genesymbol',
                'trembl',
            )

        if not uniprots and mouse:

            mouse_gs = '-'.join(
                x.capitalize()
                for x in r['official gene symbol'].split('-')
            )

            uniprots = mapping.map_name(
                mouse_gs,
                'genesymbol',
                'uniprot',
                ncbi_tax_id = 10090,
            )

            if not uniprots:

                uniprots = mapping.map_name(
                    mouse_gs,
                    'genesymbol',
                    'trembl',
                    ncbi_tax_id = 10090,
                )

            this_ncbi_tax_id = 10090

        if not uniprots and r['gene type'] == 'non-coding RNA':

            uniprots = (r['official gene symbol'],)
            entity_type = 'lncrna'

        for u in uniprots:

            result[u].add(
                record(
                    canonical_marker = r['canonical marker'] == '1',
                    cell_type = r['cell type'],
                    organ = r['organ'],
                    germ_layer = r['germ layer'],
                    entity_type = entity_type, # TODO: handle all entity types
                    human = human,
                    mouse = mouse,
                    human_sensitivity = to_float(r['sensitivity_human']),
                    mouse_sensitivity = to_float(r['sensitivity_mouse']),
                    human_specificity = to_float(r['specificity_human']),
                    mouse_specificity = to_float(r['specificity_mouse']),
                    ubiquitiousness = to_float(r['ubiquitousness index']),
                    ncbi_tax_id = this_ncbi_tax_id,
                )
            )

    return result
