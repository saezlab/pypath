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
import itertools

import pypath.inputs.common as inputs_common
import pypath.resources.urls as urls
import pypath.utils.mapping as mapping
import pypath.share.common as common
import pypath.inputs.cell as cell_input


def wojtowicz2020_raw():
    """
    Returns Supplementary Table S4 from 10.1016/j.cell.2020.07.025
    (Wojtowicz et al. 2020) as a list of tuples.
    """

    path = cell_input.cell_supplementary(
        supp_url = urls.urls['wojtowicz2020']['url'],
        article_url = urls.urls['wojtowicz2020']['article'],
    )

    content = inputs_common.read_xls(path)

    fields = content.pop(0)
    fields = [re.sub('[- ]', '_', f.lower()) for f in fields]

    Wojtowicz2020RawRecord = collections.namedtuple(
        'Wojtowicz2020RawRecord',
        fields
    )

    return [
        Wojtowicz2020RawRecord(
            *(
                float(f)
                    if 5 < i < 17 else
                f
                for i, f in enumerate(line)
            )
        )
        for line in content
    ]


def _id_translate(name):

    return mapping.map_name(name, 'genesymbol', 'uniprot')


def wojtowicz2020_interactions():

    Wojtowicz2020Interaction = collections.namedtuple(
        'Wojtowicz2020Interaction',
        ['id_a', 'id_b'],
    )

    result = []

    for rec in wojtowicz2020_raw():

        preys = _id_translate(rec.prey_gene_name)
        baits = _id_translate(rec.bait_gene_name)

        for id_a, id_b in itertools.product(preys, baits):

            result.append(Wojtowicz2020Interaction(id_a, id_b))

    return result