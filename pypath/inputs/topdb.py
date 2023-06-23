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


def topdb_annotations(ncbi_tax_id = 9606):

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
    )

    parser = etree.iterparse(c.fileobj, events = ('start', 'end'))

    result = collections.defaultdict(set)
    root = next(parser)
    used_elements = []

    for ev, elem in parser:

        if ev == 'end' and elem.tag == 'TOPDB':

            used_elements.append(elem)

            organism = elem.find('Organism').text
            organism = taxonomy.ensure_ncbi_tax_id(organism)

            if not organism:

                continue

            tag_uniprots = elem.find('./CrossRef/UniProt')

            if tag_uniprots is None:
                continue

            uniprots = [u.text for u in tag_uniprots.findall('AC')]
            uniprots = set(
                mapping.map_name0(
                    u,
                    'uniprot',
                    'uniprot',
                    ncbi_tax_id = ncbi_tax_id,
                )
                for u in uniprots
            )

            if not uniprots:
                continue

            membranes = set(
                mem
                for tag_mem in elem.findall('Membrane')
                for mem in tag_mem.text.split(';')
            )

            ntm = 0
            score = 0
            topologies = ()
            tag_topo = elem.find('Topology')

            if tag_topo is not None:
                ntm = int(tag_topo.find('Numtm').attrib['Count'])
                score = int(tag_topo.find('Reliability').text)

                topologies = set(
                    tag_reg.attrib['Loc']
                    for tag_reg in tag_topo.findall('./Regions/Region')
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
