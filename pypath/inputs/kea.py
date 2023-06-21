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
import itertools

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.utils.mapping as mapping


_resources = {
    'phosphosite': 'PhosphoSite',
    'phosphopoint': 'PhosphoPoint',
    'hprd': 'HPRD',
    'networkin': 'NetworKIN',
    'kinexus': 'Kinexus',
    'phosphoelm': 'phosphoELM',
}


def kea_interactions():


    KeaRecord = collections.namedtuple(
        'KeaRecord',
        [
            'enzyme',
            'substrate',
            'residue_type',
            'residue_offset',
            'pmid',
            'resource',
        ]
    )

    resub = re.compile(r'(\w+)_([A-Z])([0-9]+)')

    url = urls.urls['kea']['kinase_substrate']

    c = curl.Curl(url, silent = False, large = True)

    result = []

    for rec in c.result:

        rec = rec.strip().split('\t')

        site = resub.match(rec[1].strip())

        if not site:

            continue

        target, resaa, resnum = site.groups()

        e_uniprots = mapping.map_name(rec[0], 'genesymbol', 'uniprot')
        s_uniprots = mapping.map_name(target, 'genesymbol', 'uniprot')

        for enz, sub in itertools.product(e_uniprots, s_uniprots):

            result.append(
                KeaRecord(
                    enzyme = enz,
                    substrate = sub,
                    residue_type = resaa,
                    residue_offset = int(resnum),
                    pmid = rec[2].strip(),
                    resource = _resources[rec[3].strip()]
                )
            )

    return result


def kea_enzyme_substrate():

    return [
        {
            'start': None,
            'end': None,
            'instance': None,
            'substrate': rec.substrate,
            'kinase': rec.enzyme,
            'resaa': rec.residue_type,
            'resnum': rec.residue_offset,
            'references': {rec.pmid},
            'typ': 'phosphorylation',
        }
        for rec in kea_interactions()
    ]
