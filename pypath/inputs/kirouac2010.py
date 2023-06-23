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

import os
import re
import itertools
import collections

from typing import List

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.inputs.embopress as embo


def kirouac2010_interactions() -> List[tuple]:
    """
    Ligand-receptor pairs from Kirouac et al. 2010
    (https://www.embopress.org/doi/10.1038/msb.2010.71).
    """

    Kiruac2010Interaction = collections.namedtuple(
        'Kiruac2010Interaction',
        (
            'ligand',
            'receptor',
        ),
    )

    rename = re.compile(r'[A-Z]{2}[A-Z0-9][-A-Z0-9]*')
    rerange = re.compile(r'([0-9])-([0-9])')
    reslash = re.compile(r'.*?([A-Z0-9]{1,3}/[/A-Z0-9]+)')


    def get_names(s):

        names = set()
        prev = None

        for n in s.split():

            m = rename.findall(n)

            if m:
                prev = m
                m = reslash.match(n)

                if m:
                    for post in m.groups()[0].split('/'):
                        for pre in prev:
                            names.add('%s%s' % (pre, post))

                else:
                    m = rerange.match(n)

                    if m:
                        intv = m.groups()

                        for post in range(int(intv[0]), int(intv[1]) + 1):
                            for pre in prev:
                                names.add('%s%u' % (pre, post))

                    else:
                        names.update(prev)

            prev = None

        return names


    tbl = embo.embopress_supplementary(
        url = urls.urls['kirouac2010']['url'],
        init_url = urls.urls['kirouac2010']['init_url'],
        sheet = 'S12',
    )

    result = []

    for r in tbl[2:]:
        namesA = get_names(r[0])
        namesB = get_names(r[1])

        result.extend([
            Kiruac2010Interaction(*lig_rec)
            for lig_rec in itertools.product(namesA, namesB)
        ])

    return result
