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

try:
    from cStringIO import StringIO
except:
    try:
        from StringIO import StringIO
        from StringIO import StringIO as BytesIO
    except:
        from io import BytesIO
        from io import StringIO

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.inputs.common as inputs_common
import pypath.inputs.ontology as ontology


def get_switches_elm():
    """
    switches.elm is a resource containing functional switches in molecular
    regulation, in domain-motif level resolution, classified into categories
    according to their mechanism.
    """

    residue = re.compile(r'(^[A-Z])([0-9]+)')
    url = urls.urls['switches.elm']['url']
    c = curl.Curl(url, silent = False)
    data = c.result

    if data is None:

        return None

    buff = StringIO()
    buff.write(data)
    cols = {
        'intramol': 3,
        'bindingsite_a': 5,
        'bs_a_start': 6,
        'bs_a_end': 7,
        'uniprot_a': 4,
        'uniprot_b': 8,
        'bindingsite_b': 9,
        'bs_b_start': 10,
        'bs_b_end': 11,
        'affected': 12,
        'type': 13,
        'subtype': 14,
        'mechanism': 15,
        'reversible': 16,
        'outcome': 17,
        'outcomedir': 18,
        'modification': 19,
        'modsites': 20,
        'modifiers': 21,
        'effectors': 22,
        'references': 26,
    }
    subf = {
        4: 'UNIPROT:',
        8: 'UNIPROT:',
        25: ';',
        26: ';',
    }
    table = inputs_common.read_table(
        cols = cols,
        fileObject = buff,
        sep2 = subf,
        hdr = 1,
    )
    mod_ont = ontology.ontology('MOD')

    for l in table:

        if l['modification'].startswith('MOD'):

            if l['modification'] in mod_ont:

                l['modification'] = mod_ont[l['modification']]

        l['references'] = [
            x.replace('PMID:', '').strip() for x in l['references']
        ]
        l['modsites'] = [
            (m.group(2), m.group(1))
            for m in
            (
                residue.match(s.strip())
                for s in l['modsites'].split(';')
                if s
            )
        ]
        l['intramol'] = True if l['intramol'].strip() == 'TRUE' else False
        l['bs_a_start'] = [x.split(';') for x in l['bs_a_start'].strip()]
        l['bs_b_start'] = [x.split(';') for x in l['bs_b_start'].strip()]
        l['bs_a_end'] = [x.split(';') for x in l['bs_a_end'].strip()]
        l['bs_b_end'] = [x.split(';') for x in l['bs_b_end'].strip()]
        l['bindingsite_a'] = [
            x.strip()
            for x in l['bindingsite_a'].split(';')
        ]
        l['bindingsite_b'] = [
            x.strip()
            for x in l['bindingsite_b'].split(';')
        ]
        l['modifiers'] = [
            x.split(':') for x in l['modifiers'].strip().split(';')
        ]
        bs_a_ids = {}
        bs_b_ids = {}
        mod_ids = {}

        for bs in l['bindingsite_a']:

            if ':' in bs:

                bs = bs.split(':')

                if bs[0].lower() not in bs_a_ids:
                    bs_a_ids[bs[0].lower()] = []

                bs_a_ids[bs[0].lower()].append(bs[1])

        for bs in l['bindingsite_b']:

            if ':' in bs:

                bs = bs.split(':')

                if bs[0].lower() not in bs_b_ids:

                    bs_b_ids[bs[0].lower()] = []

                bs_b_ids[bs[0].lower()].append(bs[1])

        for mod in l['modifiers']:

            if ':' in mod:

                mod = mod.split(':')

                if mod[0].lower() not in mod_ids:

                    mod_ids[mod[0].lower()] = []

                mod_ids[mod[0].lower()].append(mod[1])

        l['bindingsite_a'] = bs_a_ids
        l['bindingsite_b'] = bs_b_ids
        l['modifiers'] = mod_ids

    return table
