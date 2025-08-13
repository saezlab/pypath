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

from past.builtins import xrange, range

import collections
import itertools

from lxml import etree

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.utils.mapping as mapping
import pypath.utils.taxonomy as taxonomy


def topdb_annotations(organism = 9606):

    TopdbAnnotation = collections.namedtuple(
        'TopdbAnnotation',
        ['membrane', 'topology', 'score', 'tmregions'],
    )

    result = collections.defaultdict(set)

    url = urls.urls['topdb']['url']
    c = curl.Curl(
        url,
        large = True,
        default_mode = 'rb',
        silent = False,
        slow = True,
    )

    parser = etree.iterparse(c.fileobj, events = ('start', 'end'))
    xmlns = '{https://topdb.unitmp.org}'

    result = collections.defaultdict(set)
    root = next(parser)
    used_elements = []

    for ev, elem in parser:

        if ev == 'end' and elem.tag == f'{xmlns}TOPDB':

            used_elements.append(elem)

            ncbi_tax_id = taxonomy.ensure_ncbi_tax_id(organism)

            uniprot_ac = elem.attrib['ID']

            uniprots = mapping.map_name(
                uniprot_ac,
                'uniprot-entry',
                'uniprot',
                ncbi_tax_id = ncbi_tax_id,
            )

            if not uniprots:
                continue

            tag_mems = elem.find(f'{xmlns}Membranes')
            membranes = {'Unknown'}

            if tag_mems is not None:

                membranes = set(
                    mem
                    for tag_mem in tag_mems.findall(f'{xmlns}Membrane')
                    for mem in tag_mem.text.split(';')
                )

            ntm = 0
            score = 0
            topologies = ()
            tag_topo = elem.find(f'{xmlns}Topology')

            if tag_topo is not None:
                ntm = int(tag_topo.find(f'{xmlns}Numtm').attrib['Count'])
                score = float(tag_topo.find(f'{xmlns}Reliability').text)

                topologies = set(
                    tag_reg.attrib['Loc']
                    for tag_reg in tag_topo.findall(f'./{xmlns}Regions/{xmlns}Region')
                )

            if not membranes:
                membranes = (None,)

            if not topologies:
                topologies = (None,)

            for topology, membrane, uniprot in itertools.product(
                topologies,
                membranes,
                uniprots,
            ):

                if uniprot is None:

                    continue

                result[uniprot].add(
                    TopdbAnnotation(
                        membrane = membrane,
                        topology = topology,
                        tmregions = ntm,
                        score = score,
                    )
                )

        # removing used elements to keep memory low
        if len(used_elements) > 2000:
            for _ in xrange(1000):
                e = used_elements.pop(0)
                e.clear()

    # closing the XML
    c.fileobj.close()
    del c

    return dict(result)
