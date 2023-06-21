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

import collections

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.utils.mapping as mapping
import pypath.inputs.common as inputs_common


def matrisome_annotations(organism = 9606):
    """
    Downloads MatrisomeDB 2.0, a database of extracellular matrix proteins.
    Returns dict where keys are UniProt IDs and values are tuples of
    classes, subclasses and notes.
    """

    MatrisomeAnnotation = collections.namedtuple(
        'MatrisomeAnnotation',
        ['mainclass', 'subclass', 'subsubclass']
    )


    tax_names = {
        10090: ('Murine', 'mm'),
        9606:  ('Human',  'hs'),
    }

    url = urls.urls['matrisome']['url_rescued'] % tax_names[organism][1]
    c = curl.Curl(url, large = True, silent = False)
    xlsname = c.fname
    del(c)
    raw = inputs_common.read_xls(xlsname)[1:]

    result = collections.defaultdict(set)

    for r in raw:

        uniprots = set(r[7].split(':'))
        uniprots.discard('')

        if not uniprots:
            continue

        uniprots = mapping.map_names(uniprots, 'uniprot', 'uniprot')

        for uniprot in uniprots:

            result[uniprot].add(
                MatrisomeAnnotation(
                    mainclass = r[0].strip(),
                    subclass = r[1].strip(),
                    subsubclass = r[10].strip() or None,
                )
            )

    return dict(result)


def __matrisome_annotations_2():
    """
    This I made only to find out why certain proteins are missing from this
    output. I will contact Matrisome people to ask why.
    """

    url = urls.urls['matrisome']['url_dl']
    c = curl.Curl(url, large = True, silent = False)

    _ = next(c.result)

    return set(r.split(',')[1] for r in c.result)
