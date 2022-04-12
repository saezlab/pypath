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

import pypath.share.common as common


def mitab_field_list(field):

    return common.uniq_list(
        map(
            lambda x: x.split('(')[1][:-1],
            field.split('|')
        )
    )


def mitab_field_uniprot(field):

    uniprots = list(
        filter(
            lambda x: len(x) == 2 and x[0] == 'uniprotkb',
            map(
                lambda x: x.split(':'),
                field.split('|')
            )
        )
    )

    return uniprots[0][1] if uniprots else None
