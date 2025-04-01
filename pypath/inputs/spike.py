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

import itertools
import collections

from typing import Dict, List

import xml.etree.cElementTree as ET

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.utils.mapping as mapping
import pypath.internals.intera as intera


def spike_interactions(min_confidence: int = 2) -> List[tuple]:
    """
    Args
        min_confidence:
            Confidence (integrity) levels in SPIKE span from 1 to 4, 1
            being the highest confidence.
    """

    url = urls.urls['spike']['url']
    c = curl.Curl(
        url,
        silent = False,
        large = True,
        files_needed = ['LatestSpikeDB.xml'],
    )
    spikexml = c.result

    xml = ET.parse(spikexml['LatestSpikeDB.xml'])

    xmlroot = xml.getroot()

    # iterating genes

    bblock = xmlroot.find('BuildingBlock')
    rblock = xmlroot.find('RegulationBlock')
    iblock = xmlroot.find('InteractionBlock')

    genes = {}
    result = []

    SpikeInteraction = collections.namedtuple(
        'SpikeInteraction',
        (
            'entrez_a',
            'genesymbol_a',
            'entrez_b',
            'genesymbol_b',
            'directed',
            'pmids',
            'integrity',
            'effect',
            'assay',
            'data_source',
            'description',
            'mechanism',
            'regulation',
        )
    )

    Gene = collections.namedtuple('Gene', ('entrez', 'genesymbol', 'type'))

    for gene in bblock.findall('Gene'):

        sy = '' if 'name' not in gene.attrib else gene.attrib['name']
        genes[gene.attrib['id']] = [
            Gene(
                entrez = gene.find('XRef').attrib['id'],
                genesymbol = sy,
                type = 'gene',
            )
        ]

    for grp in bblock.findall('Group'):

        if grp.attrib['type'] == 'Complex':

            members = [m.attrib['ref'] for m in grp.findall('Member')]

            try:

                if all(
                    m in genes and
                    genes[m][0].type != 'complex'
                    for m in members
                ):

                    uniprots = [
                        mapping.map_name(genes[m][0].entrez, 'entrez', 'uniprot')
                        for m in members
                    ]

                    genes[grp.attrib['id']] = [
                        Gene(
                            entrez = cplex,
                            genesymbol = cplex,
                            type = 'complex',
                        )
                        for cplex in
                        (
                            intera.Complex(
                                name = grp.attrib['name'],
                                components = ups,
                                sources = 'SPIKE',
                            )
                            for ups in itertools.product(*uniprots)
                        )
                    ]

            except Exception as e:

                for m in members:
                    # print(m)
                    # print(genes[m])
                    pass
                raise e

    for i in itertools.chain(
            rblock.findall('Regulation'),
            iblock.findall('Interaction'),
        ):

        regulation = i.tag == 'Regulation'
        src_tag = 'Source' if regulation else 'ProteinA'
        tgt_tag = 'PhysicalTarget' if regulation else 'ProteinB'

        ds = i.attrib['dataSource']
        itg = i.attrib['integrity']
        eff = i.attrib.get('effect', '')
        mec = i.attrib.get('mechanism', '')
        src = i.find(src_tag).attrib['ref']
        tgt = i.find(tgt_tag).attrib['ref']
        dcd = str(int(regulation))
        dsc = (
            ''
                if i.find('Description') is None else
            i.find('Description').text.replace('\n', ' ')
        )
        asy = (
            ''
                if 'biologicalAssay' not in i.attrib else
            i.attrib['biologicalAssay']
        )
        refs = i.findall('Reference')
        pmids = [r.attrib['pmid'] for r in refs]

        if src in genes and tgt in genes:

            if int(itg) <= min_confidence:

                for _src, _tgt in itertools.product(genes[src], genes[tgt]):

                    result.append(
                        SpikeInteraction(
                            entrez_a = _src.entrez,
                            genesymbol_a = _src.genesymbol,
                            entrez_b = _tgt.entrez,
                            genesymbol_b = _tgt.genesymbol,
                            directed = dcd,
                            pmids = ';'.join(pmids),
                            integrity = itg,
                            effect = eff,
                            assay = asy,
                            data_source = ds,
                            description = dsc,
                            mechanism = mec,
                            regulation = regulation,
                        )
                    )

    return result


def spike_complexes(min_confidence: int = 2) -> Dict[str, intera.Complex]:

    interactions = spike_interactions(min_confidence = min_confidence)

    complexes = [
        getattr(i, attr)
        for i in interactions
        for attr in ('entrez_a', 'entrez_b')
        if isinstance(getattr(i, attr), intera.Complex)
    ]

    return dict((cplx.__str__(), cplx) for cplx in complexes)
