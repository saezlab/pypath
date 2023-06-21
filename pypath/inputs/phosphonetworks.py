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

import pypath.share.curl as curl
import pypath.resources.urls as urls


def phosphonetworks_enzyme_substrate():

    result = []
    reres = re.compile(r'([A-Z])([0-9]+)')
    non_digit = re.compile(r'[^\d.-]+')
    motre = re.compile(r'(-*)([A-Za-z]+)(-*)')
    url = urls.urls['phosnw']['url']
    c = curl.Curl(url, silent = False)
    data = c.result

    if data is None:

        return

    data = data.split('\n')

    for l in data:

        if l.startswith('>'):

            substrate = l[1:].strip()

        elif len(l.split('\t')) >= 4:

            l = [x.strip() for x in l.split('\t')]
            res = reres.match(l[1])
            resnum = int(non_digit.sub('', res.groups()[1]))
            mot = motre.match(l[0])

            if mot:
                start = resnum - 7 + len(mot.groups()[0])
                end = resnum + 7 - len(mot.groups()[2])
                instance = l[0].replace('-', '').upper()
            else:
                start = None
                end = None
                instance = l[0]
            result.append({
                'instance': instance,
                'kinase': l[2],
                'resaa': res.groups()[0],
                'resnum': resnum,
                'score': float(non_digit.sub('', l[3])),
                'substrate': substrate,
                'start': start,
                'end': end
            })

    return result


def phosphonetworks_interactions():

    result = []
    data = phosphonetworks_enzyme_substrate()
    for l in data:
        result.append((l['kinase'], l['substrate']))

    return [list(x) for x in list(set(result))]
