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
import itertools

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.utils.mapping as mapping
import pypath.utils.taxonomy as taxonomy
import pypath.share.common as common


def depod_interactions(organism = 9606):

    url = urls.urls['depod']['urls'][1]
    c = curl.Curl(url, silent = False, large = True, encoding = 'iso-8859-1')
    data = c.result
    result = []
    i = []
    lnum = 0

    for l in data:

        if lnum == 0:
            lnum += 1
            continue
        l = l.replace('\n', '').replace('\r', '')
        l = l.split('\t')
        specA = int(l[9].split(':')[1].split('(')[0])
        specB = int(l[10].split(':')[1].split('(')[0])

        if organism is None or (specA == organism and specB == organism):

            pm = l[8].replace('pubmed:', '')
            sc = l[14].replace('curator score:', '')
            ty = l[11].split('(')[1].replace(')', '')
            l = [l[0], l[1]]
            interaction = ()

            for ll in l:

                ll = ll.split('|')
                uniprot = ''
                for lll in ll:
                    nm = lll.split(':')
                    u = nm[1].strip()
                    if nm[0] == 'uniprotkb' and len(u) == 6:
                        uniprot = u
                interaction += (uniprot, )

            interaction += (pm, sc, ty)
            if len(interaction[0]) > 1 and len(interaction[1]) > 1:
                i.append(interaction)

        lnum += 1

    return i


def depod_enzyme_substrate(organism = 9606):

    result = []

    reunip = re.compile(r'uniprotkb:([A-Z0-9]+)')
    reptm = re.compile(r'([A-Z][a-z]{2})-([0-9]+)')
    repmidsep = re.compile(r'[,|]\s?')

    url = urls.urls['depod']['urls'][0]
    c = curl.Curl(url, silent = False, encoding = 'ascii')
    data = c.result
    data = [x.split('\t') for x in data.split('\n')]
    del data[0]

    url_mitab = urls.urls['depod']['urls'][1]
    c_mitab = curl.Curl(url_mitab, silent = False, encoding = 'iso-8859-1')
    data_mitab = c_mitab.result
    data_mitab = [x.split('\t') for x in data_mitab.split('\n')]
    del data_mitab[0]

    for i, l in enumerate(data):

        if (
            len(l) > 6 and
            l[2] == 'protein substrate' and
            taxonomy.ensure_ncbi_tax_id(
                l[3].split('(')[0].strip()
            ) == organism and
            l[4].strip() != 'N/A'
        ):

            enzyme_uniprot = reunip.search(data_mitab[i][0]).groups()[0]
            substrate_uniprot = reunip.search(data_mitab[i][1]).groups()[0]

            for enzyme_up, substrate_up in itertools.product(
                    mapping.map_name(
                        enzyme_uniprot,
                        'uniprot',
                        'uniprot'
                    ),
                    mapping.map_name(
                        substrate_uniprot,
                        'uniprot',
                        'uniprot'
                    ),
                ):

                for resaa, resnum in reptm.findall(l[4]):

                    resnum = int(resnum)
                    resaa = (
                        common.aminoa_3_to_1_letter[resaa]
                            if resaa in common.aminoa_3_to_1_letter else
                        resaa
                    )

                    result.append({
                        'instance': None,
                        'kinase': enzyme_up,
                        'resaa': resaa,
                        'resnum': resnum,
                        'references': repmidsep.split(l[6].strip()),
                        'substrate': substrate_up,
                        'start': None,
                        'end': None,
                        'typ': 'dephosphorylation',
                    })

    return result
