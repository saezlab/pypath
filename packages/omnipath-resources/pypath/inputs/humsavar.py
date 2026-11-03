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
import collections

import pypath.resources.urls as urls
import pypath.share.curl as curl

redesc = re.compile(
    r'([^\(^\[]+[^\s^\(^\[])\s?'
    r'(?:\((\w+)\))?\s?'
    r'(?:\[(MIM:\d+)\])?'
)
reempty = re.compile(r'^-?\s?$')


def _parse_desc(desc):

        return (
            (None,) * 3
                if reempty.match(desc) else
            redesc.match(desc).groups()
        )


def uniprot_variants() -> dict[str, set[tuple]]:
    """
    Retrieves all human missense variants annotated in UniProtKB/Swiss-Prot.

    Returns:
        Drug attributes in below as a list of named tuples.
    """

    UniprotVariant = collections.namedtuple(
            'UniprotVariant',
            (
                'genesymbol',
                'ftid',
                'aa_change',
                'variant_category',
                'dbsnp',
                'disease',
                'disease_symbol',
                'disease_omim'
            ),
        )

    url = urls.urls['humsavar']['url']
    c = curl.Curl(url, large=True, silent=False)

    result = collections.defaultdict(set)
    respace = re.compile(r'\s+')

    # skipping data description information in txt file
    for r in c.result:

        if r.startswith('_'):

            break

    for line in c.result:

        line = respace.split(line.strip())

        if len(line) == 1:

            break

        disease, symbol, omim = _parse_desc(' '.join(line[6:]))

        variant = UniprotVariant(
                line[0],
                *line[2:6],
                disease = disease,
                disease_symbol = symbol,
                disease_omim = omim,
        )

        result[line[1]].add(variant)

    return dict(result)
