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

from future.utils import iteritems

import re
import bs4

import pypath.utils.taxonomy as taxonomy
import pypath.share.curl as curl
import pypath.resources.urls as urls


def dbptm_enzyme_substrate(organism = 9606):
    """
    Downloads enzyme-substrate interactions from dbPTM.
    Returns list of dicts.
    """

    if organism is None:
        _organism = None
    elif organism in taxonomy.dbptm_taxids:
        _organism = taxonomy.dbptm_taxids[organism]
    else:
        sys.stdout.write('\t:: Unknown organism: `%u`.\n' % organism)
        return []

    url = urls.urls['dbptm']['old_table']
    c = curl.Curl(url, silent = False, large = True)
    data = []

    hdr = next(c.result).strip().split('\t')

    for l in c.result:

        l = l.strip().split('\t')

        data.append(dict(
            (
                key,
                (
                    None
                        if val == '' else
                    val.split(';')
                        if key in {'references', 'kinase'} else
                    int(val)
                        if val.isdigit() else
                    val
                )
            )
            for key, val in zip(hdr, l)
        ))

    return data


def dbptm_enzyme_substrate_old(organism = 9606):
    """
    Downloads enzyme-substrate interactions from dbPTM.
    Returns list of dicts.
    """

    if organism is None:
        _organism = None
    elif organism in taxonomy.dbptm_taxids:
        _organism = taxonomy.dbptm_taxids[organism]
    else:
        sys.stdout.write('\t:: Unknown organism: `%u`.\n' % organism)
        return []

    result = []
    byre = re.compile(r'.*by\s([A-Za-z0-9\s]+)\.*')
    andre = re.compile(r',|and')
    non_digit = re.compile(r'[^\d.-]+')

    for url in urls.urls['dbptm']['urls']:

        c = curl.Curl(url, silent = False)
        extra = c.result

        for k, data in iteritems(extra):

            data = [x.split('\t') for x in data.split('\n')]

            for l in data:

                if len(l) > 8:

                    if _organism:
                        mnemonic = l[0].split('_')[1].strip()
                        if mnemonic != _organism:
                            continue

                    resnum = int(non_digit.sub('', l[2]))

                    ptm = ({
                        'substrate': l[1],
                        'typ': l[7].lower(),
                        'resaa': l[8][6],
                        'resnum': resnum,
                        'instance': l[8].strip(),
                        'references': l[4].split(';'),
                        'databases': (l[5].split()[0],),
                        'kinase': None if byre.match(l[3]) is None else [
                            i.strip()
                            for i in andre.split(
                                byre.match(l[3]).groups(1)[0])
                        ],
                        'start': resnum - 6,
                        'end': resnum + 6,
                    })

                    if ptm['kinase'] is not None:

                        if 'autocatalysis' in ptm['kinase']:

                            ptm['kinase'].append(ptm['substrate'])
                            ptm['kinase'].remove('autocatalysis')

                        ptm['kinase'] = [
                            k.replace('host', '').strip()
                            for k in ptm['kinase']
                        ]

                        ptm['kinase'] = [
                            k for k in ptm['kinase'] if len(k) > 0
                        ]

                        if len(ptm['kinase']) == 0:
                            ptm['kinase'] = None

                    result.append(ptm)

    return result


def dbptm_interactions():

    result = []
    data = dbptm_enzyme_substrate()
    for r in data:

        if r['kinase'] is not None:

            for src in r['kinase']:

                result.append([
                    src,
                    r['substrate'],
                    ';'.join(
                        i
                        for i in r['references']
                        if i != '-'
                    ),
                ])

    return result
