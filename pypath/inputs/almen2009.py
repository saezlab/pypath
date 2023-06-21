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

import pypath.inputs.common as inputs_common
import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.utils.mapping as mapping


def almen2009_annotations():


    resep = re.compile(r'[;/]')


    Almen2009Annotation = collections.namedtuple(
        'Almen2009Annotation',
        [
            'mainclass',
            'classes',
            'phobius_secreted',
            'phobius_transmembrane',
            'sosui_transmembrane',
            'tmhmm_transmembrane',
        ]
    )


    url = urls.urls['almen2009']['url']

    c = curl.Curl(url, silent = False, large = True)

    xls = c.fileobj
    xlsfile = xls.name
    xls.close()
    tbl = inputs_common.read_xls(xlsfile, sheet = 'Data')[1:]

    result = collections.defaultdict(set)

    for row in tbl:

        uniprots = mapping.map_name(row[0], 'ipi', 'uniprot')

        mainclass = row[2]
        classes = row[3].replace('KInase', 'Kinase')
        classes = tuple(sorted(resep.split(classes)))
        phobius_transmembrane = int(float(row[5]))
        phobius_secreted = row[6] == 'Y'
        sosui_transmembrane = int(float(row[8])) if row[8] != 'ERROR' else 0
        tmhmm_transmembrane = int(float(row[10]))

        for uniprot in uniprots:

            result[uniprot].add(
                Almen2009Annotation(
                    mainclass = mainclass,
                    classes = classes,
                    phobius_secreted = phobius_secreted,
                    phobius_transmembrane = phobius_transmembrane,
                    sosui_transmembrane = sosui_transmembrane,
                    tmhmm_transmembrane = tmhmm_transmembrane,
                )
            )

    return dict(result)
