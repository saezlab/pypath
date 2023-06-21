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


def get_hprd(in_vivo = True):
    """
    Downloads and preprocesses HPRD data.
    """

    url = urls.urls['hprd_all']['url_rescued']
    files = [urls.urls['hprd_all']['ptm_file']]
    c = curl.Curl(url, silent = False, files_needed = files)
    data = c.result

    if len(data) == 0:
        return []

    data = [l.split('\t') for l in data[files[0]].split('\n')][:-1]

    if in_vivo:
        data = [i for i in data if 'in vivo' in i[9].split(';')]

    return data


def hprd_interactions(in_vivo = True):
    """
    Processes HPRD data and extracts interactions.
    Returns list of interactions.
    """
    return [i for i in get_hprd(in_vivo = in_vivo) if i[6] != '-']


def hprd_interactions_htp():

    url = urls.urls['hprd_all']['url_rescued']
    fname = urls.urls['hprd_all']['int_file']
    c = curl.Curl(url, silent = False, large = True, files_needed = [fname])

    return list(
        map(
            lambda l: l.split('\t'),
            c.result[fname].read().decode('ascii').split('\n')
        )
    )


def hprd_enzyme_substrate(in_vivo = True):
    """
    Processes HPRD data and extracts PTMs.
    Returns list of kinase-substrate interactions.
    """

    ptms = []
    non_digit = re.compile(r'[^\d]+')
    data = get_hprd(in_vivo = in_vivo)

    for ptm in data:

        if ptm[6] != '-':
            resnums = [
                int(nn)
                for nn in [non_digit.sub('', n) for n in ptm[4].split(';')]
                if len(nn) > 0
            ]

            for resnum in resnums:

                modtype = ptm[8].lower().replace('proteolytic', '').strip()

                ptms.append({
                    'resaa': ptm[5],
                    'resnum': resnum,
                    'typ': modtype,
                    'references': ptm[10].split(','),
                    'kinase': ptm[6],
                    'substrate_refseqp': ptm[3],
                    'substrate': ptm[1],
                    'start': max(resnum - 7, 1),
                    'end': resnum + 7,
                    'instance': None
                })

    return ptms
