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

import collections
import json

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.utils.mapping as mapping
import pypath.share.progress as progress


def membranome_annotations():

    membr_url = urls.urls['membranome']['baseurl'] % ('membranes', '')
    c = curl.Curl(membr_url, large = True, silent = False)
    membr_data = json.loads(c.fileobj.read())
    del c

    membr = dict((m['id'], m) for m in membr_data['objects'])

    page = 1
    prot_all = []

    prg = progress.Progress(7, 'Downloading Membranome', 1)

    while True:

        prg.step()

        prot_url = urls.urls['membranome']['baseurl'] % (
            'proteins',
            '?pageSize=1000&pageNum=%u' % page,
        )
        c = curl.Curl(prot_url, large = True, silent = True)
        prot = json.loads(c.fileobj.read())

        prot_all.extend(prot['objects'])

        if prot['page_end'] >= prot['total_objects']:
            break

        page = prot['page_num'] + 1

    prg.terminate()

    for p in prot_all:

        uniprots = mapping.map_name(p['uniprotcode'], 'uniprot', 'uniprot')

        for uniprot in uniprots:

            yield (
                uniprot,
                membr[p['membrane_id']]['name'],
                membr[p['membrane_id']]['topology_in']
                    if p['topology_show_in'] else
                membr[p['membrane_id']]['topology_out'],
            )
