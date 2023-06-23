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
import collections
import bs4

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.utils.mapping as mapping
import pypath.utils.reflists as reflists


def tcdb_families():

    retag = re.compile(r'<.*>')
    rethe = re.compile(r'^[tT]he (.*) (?:Super)?[Ff]amily')

    url = urls.urls['tcdb']['url_families']

    c = curl.Curl(url, large = False, silent = False)

    lines = bs4.BeautifulSoup(c.result, features = 'lxml').find('p').text

    return dict(
        (
            tcid,
            rethe.sub(r'\g<1>', family.replace('\t', ' '))
        )
        for tcid, family in
        (
            retag.sub('', line.strip()).split('\t', maxsplit = 1)
            for line in lines.strip().split('\n')
        )
    )


def tcdb_classes():

    refam = re.compile(r'(\d\.[A-Z]\.\d+)')
    retab = re.compile(r'\t+')

    url = urls.urls['tcdb']['url_acc2tc']

    c = curl.Curl(url, large = True, silent = False)

    result = {}

    for line in c.result:

        if not line:

            continue

        ac, tc = retab.split(line.rstrip())
        family = refam.search(tc).groups()[0]

        result[ac] = (tc, family)

    return result


def tcdb_annotations(organism = 9606):


    TcdbAnnotation = collections.namedtuple(
        'TcdbAnnotation',
        [
            'family',
            'tcid',
        ]
    )


    families = tcdb_families()
    classes = tcdb_classes()
    result = collections.defaultdict(set)

    for ac, (tc, family) in iteritems(classes):

        uniprots = mapping.map_name(
            ac,
            'uniprot',
            'uniprot',
            ncbi_tax_id = organism,
        )

        for uniprot in uniprots:

            if reflists.check(uniprot, 'uniprot', ncbi_tax_id = organism):

                result[uniprot].add(
                    TcdbAnnotation(
                        family = families[family],
                        tcid = tc,
                    )
                )

    return dict(result)
