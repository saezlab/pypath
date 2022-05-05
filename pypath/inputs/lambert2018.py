#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2022
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: Dénes Türei (turei.denes@gmail.com)
#           Nicolàs Palacio
#           Sebastian Lobentanzer
#           Erva Ulusoy
#           Olga Ivanova
#           Ahmet Rifaioglu
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

import re
import collections

import pypath.inputs.common as inputs_common
import pypath.resources.urls as urls
import pypath.utils.mapping as mapping
import pypath.inputs.cell as cell_input
import pypath.share.common as common


def lambert2018_s1_raw():

    path = cell_input.cell_supplementary(
        supp_url = urls.urls['lambert2018']['s1'],
        article_url = urls.urls['lambert2018']['article'],
    )

    content = inputs_common.read_xls(path, sheet = 1)

    h0, h1 = content.pop(0), content.pop(0)
    h1[3] = h0[3]

    names = ['%s_%s' % n for n in zip()]

    record = collections.namedtuple(
        'Lambert2018Raw',
        [
            nn for nn in (
                re.sub('[- ?;:]', '_', n).lower().strip('_ ')
                for n in h1
            )
            if nn
        ]
    )

    nfields = len(record._fields)

    return [
        record(*(common.try_bool(common.try_float(f)) for f in r[:nfields]))
        for r in content
    ]
