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
import csv
import collections

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.utils.mapping as mapping
import pypath.utils.taxonomy as taxonomy


def opm_annotations(organism = 9606):

    reparentheses = re.compile(r'\((.*)\)')
    regenesymbol  = re.compile(r' ([A-Z0-9]{3,}) ')


    def get_dict(name):

        result = {}
        url = urls.urls['opm'][name]
        c = curl.Curl(url, large = True, silent = False)
        data = csv.DictReader(c.result, delimiter = ',')

        for rec in data:
            result[rec['id']] = rec['name']

        return result


    OpmAnnotation = collections.namedtuple(
        'OpmAnnotation',
        ['membrane', 'family', 'transmembrane'],
    )


    result = collections.defaultdict(set)

    organism_name = (
        taxonomy.phosphoelm_taxids[organism]
            if organism in taxonomy.phosphoelm_taxids else
        None
    )

    types = get_dict('types')
    families = get_dict('families')

    url = urls.urls['opm']['proteins']
    c = curl.Curl(url, silent = False, large = True)

    data = csv.DictReader(c.result, delimiter = ',')

    for rec in data:

        if organism_name and rec['species_name_cache'] != organism_name:

            continue

        name = rec['name']

        names = [
            name,
            name.split('(')[0],
            name.split(',')[0],
        ]

        m = reparentheses.search(name)

        if m:
            names.append(m.groups()[0])

        genesymbols = regenesymbol.findall(name)

        for this_name in names:
            uniprot = mapping.map_name0(this_name, 'protein-name', 'uniprot')

            if uniprot:
                break

        if not uniprot:
            for gs in genesymbols:
                uniprot = (
                    mapping.map_name0(this_name, 'genesymbol', 'uniprot')
                )

                if uniprot:
                    break

        if not uniprot:
            continue

        result[uniprot].add(
            OpmAnnotation(
                membrane = rec['membrane_name_cache'],
                family = rec['family_name_cache'],
                transmembrane = types[rec['type_id']] == 'Transmembrane',
            )
        )

    return dict(result)
