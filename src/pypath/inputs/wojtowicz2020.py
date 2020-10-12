#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2020
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#                  Nicolàs Palacio
#                  Olga Ivanova
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

import collections

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

    return content


    #Wojtowicz2020RawRecord = collections.namedtuple(
        #'Wojtowicz2020RawRecord',
        #content[0]
    #)

    #return [
        #Wojtowicz2020RawRecord(
            #*(line[:2] + [int(float(n)) for n in line[2:]])
        #)
        #for line in
        #content[1:]
    #]