#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2020
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#                  Nicolàs Palacio
#                  Olga Ivanova
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

import xml.etree.cElementTree as ET

import pypath.share.curl as curl
import pypath.resources.urls as urls


def spike_interactions(high_confidence = True):

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

    for gene in bblock.findall('Gene'):

        sy = '' if 'name' not in gene.attrib else gene.attrib['name']
        genes[gene.attrib['id']] = (gene.find('XRef').attrib['id'], sy)

    for reg in rblock.findall('Regulation'):

        ds = reg.attrib['dataSource']
        itg = reg.attrib['integrity']
        eff = reg.attrib['effect']
        mec = reg.attrib['mechanism']
        src = reg.find('Source').attrib['ref']
        tgt = reg.find('PhysicalTarget').attrib['ref']
        dcd = "1"
        dsc = (
            ''
                if reg.find('Description') is None else
            reg.find('Description').text.replace('\n', ' ')
        )
        asy = (
            ''
                if 'biologicalAssay' not in reg.attrib else
            reg.attrib['biologicalAssay']
        )
        refs = reg.findall('Reference')
        pmids = []

        for r in refs:
            pmids.append(r.attrib['pmid'])

        if src in genes and tgt in genes:

            if itg == '1' or not high_confidence:

                result.append([
                    genes[src][0],
                    genes[src][1],
                    genes[tgt][0],
                    genes[tgt][1],
                    dcd,
                    ';'.join(pmids),
                    itg,
                    eff,
                    asy,
                    ds,
                    dsc,
                    mec,
                ])

    for ict in iblock.findall('Interaction'):

        ds = ict.attrib['dataSource']
        itg = ict.attrib['integrity']
        eff = '' if 'effect' not in ict.attrib else ict.attrib['effect']
        src = ict.find('ProteinA').attrib['ref']
        tgt = ict.find('ProteinB').attrib['ref']
        dcd = "0"
        mec = ""
        dsc = (
            ''
                if ict.find('Description') is None else
            ict.find('Description').text.replace('\n', ' ')
        )
        asy = (
            ''
                if 'biologicalAssay' not in ict.attrib else
            ict.attrib['biologicalAssay']
        )
        refs = ict.findall('Reference')
        pmids = []

        for r in refs:
            pmids.append(r.attrib['pmid'])

        if src in genes and tgt in genes:

            if itg == '1' or not high_confidence:

                result.append([
                    genes[src][0],
                    genes[src][1],
                    genes[tgt][0],
                    genes[tgt][1],
                    dcd,
                    ';'.join(pmids),
                    itg,
                    eff,
                    asy,
                    ds,
                    dsc,
                    mec,
                ])

    return result
