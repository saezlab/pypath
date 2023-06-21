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


def mimp_enzyme_substrate():

    db_names = {
        'PhosphoSitePlus': 'PhosphoSite',
        'PhosphoELM': 'phosphoELM',
    }

    result = []
    non_digit = re.compile(r'[^\d.-]+')
    motre = re.compile(r'(-*)([A-Za-z]+)(-*)')
    url = urls.urls['mimp']['url']
    c = curl.Curl(url, silent = False)
    data = c.result
    kclass = get_kinase_class()
    if data is None:
        return None
    data = [x.split('\t') for x in data.split('\n')]
    del data[0]

    for l in data:

        if len(l) > 6 and len(l[2]) > 0:

            kinases = l[2].split(';')
            kinases_gnames = []

            for k in kinases:
                if k.endswith('GROUP'):
                    grp = k.split('_')[0]
                    if grp in kclass['groups']:
                        kinases_gnames += kclass['groups'][grp]
                    elif grp in kclass['families']:
                        kinases_gnames += kclass['families'][grp]
                    elif grp in kclass['subfamilies']:
                        kinases_gnames += kclass['subfamilies'][grp]
                else:
                    kinases_gnames.append(k)

            mot = motre.match(l[4])

            for k in kinases_gnames:

                resaa = l[4][7]
                resnum = int(non_digit.sub('', l[3]))

                if mot:
                    start = resnum - 7 + len(mot.groups()[0])
                    end = resnum + 7 - len(mot.groups()[2])
                    instance = l[4].replace('-', '').upper()
                else:
                    start = None
                    end = None
                    instance = l[4]

                databases = [
                    db_names[db] if db in db_names else db
                    for db in l[6].split(';')
                ]

                result.append({
                    'instance': instance,
                    'kinase': k.upper(),
                    'resaa': resaa,
                    'resnum': resnum,
                    'npmid': int(non_digit.sub('', l[5])),
                    'substrate_refseq': l[1],
                    'substrate': l[0],
                    'start': start,
                    'end': end,
                    'databases': databases,
                })

    return result


def get_kinase_class():

    result = {'groups': {}, 'families': {}, 'subfamilies': {}, 'kinases': {}}
    tabs = re.compile(r'[\t]{3,}')
    reps = re.compile(r'ps[0-9]*$')
    url = urls.urls['kinclass']['rescued']
    c = curl.Curl(url, silent = False)
    data = c.result
    data = tabs.sub('', data)
    data = [x.split('\t') for x in data.split('\n')]
    data = data[9:]

    for l in data:

        if len(l) > 4:

            kinase = reps.sub('', l[0])
            group = l[2]
            family = l[3]
            subfamily = l[4]
            if group not in result['groups']:
                result['groups'][group] = []
            result['groups'][group].append(kinase)
            if family not in result['families']:
                result['families'][family] = []
            result['families'][family].append(kinase)
            if subfamily not in result['subfamilies']:
                result['subfamilies'][subfamily] = []
            result['subfamilies'][subfamily].append(kinase)
            result['kinases'][kinase] = {
                'group': group,
                'family': family,
                'subfamily': subfamily
            }

    return result


def mimp_interactions():

    result = []
    mimp = mimp_enzyme_substrate()

    for m in mimp:

        result.append([m['kinase'], m['substrate'], m['databases']])

    return result
