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

import xml.etree.cElementTree as ET

import pypath.share.curl as curl
import pypath.resources.urls as urls


def mppi_interactions(organism = 9606):

    url = urls.urls['mppi']['url_rescued']
    c = curl.Curl(url, silent = False, large = True)
    xmlfile = c.gzfile

    prefix = '{net:sf:psidev:mi}'

    result = []

    xml = ET.parse(xmlfile)

    xmlroot = xml.getroot()

    ilist = xmlroot[0][1]

    proteinInteractor = './/%sproteinInteractor' % prefix
    _organism = './/%sorganism' % prefix
    organism = '%u' % organism
    ncbiTaxId = 'ncbiTaxId'
    primaryRef = './/%sprimaryRef' % prefix
    bibref = './/%sbibref' % prefix
    interactionDetection = './/%sinteractionDetection' % prefix
    shortLabel = './/%sshortLabel' % prefix
    fullName = './/%sfullName' % prefix

    for i in ilist:

        _proteins = i.findall(proteinInteractor)

        if (
            len(_proteins) == 2 and
            (
                organism is None or
                (
                    _proteins[0].findall(
                        _organism
                    )[0].attrib[ncbiTaxId] == organism and
                    _proteins[1].findall(
                        _organism
                    )[0].attrib[ncbiTaxId] == organism
                )
            )
        ):

            pmids = []
            pms = i.findall(bibref)[0].findall(primaryRef)

            for pm in pms:

                if 'id' in pm.attrib:

                    pmids.append(pm.attrib['id'])

            meths = []
            dets = i.findall(interactionDetection)[0].findall(shortLabel)

            for m in dets:

                meths.append(m.text)

            proteins = []

            for prot in _proteins:

                thisP = {}

                if 'id' in prot.findall(primaryRef)[0].attrib:

                    thisP['u'] = prot.findall(primaryRef)[0].attrib['id']

                else:

                    thisP['u'] = ''

                thisP['nt'] = prot.findall(primaryRef)[0].attrib['db']
                thisP['gn'] = prot.findall(fullName)[0].text
                thisP['o'] = prot.findall(_organism)[0].attrib[ncbiTaxId]
                proteins.append(thisP)

            result.append([
                ';'.join(pmids),
                ';'.join(pmids),
                proteins[0]['u'],
                proteins[0]['nt'],
                proteins[0]['gn'],
                proteins[0]['o'],
                proteins[1]['u'],
                proteins[1]['nt'],
                proteins[1]['gn'],
                proteins[1]['o'],
            ])

    return result
