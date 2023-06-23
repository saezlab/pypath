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

from past.builtins import xrange, range

import csv
import collections

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.common as common
import pypath.share.settings as settings
import pypath.utils.mapping as mapping


def intogen_annotations():
    """
    Returns a list of cancer driver genes with their annotations,
    according to the IntOGen database.
    """

    IntogenAnnotation = collections.namedtuple(
        'IntogenAnnotation',
        [
            'type',
            'role',
            'curated',
            'oncodrive_role_prob',
        ],
    )


    url = urls.urls['intogen']['db2014_2']

    with settings.context(curl_connect_timeout = 100):

        c = curl.Curl(
            url,
            large = True,
            silent = False,
            files_needed = ['Drivers_type_role.tsv'],
            compr = 'zip',
        )

    for _ in xrange(7):

        __ = c.result['Drivers_type_role.tsv'].readline()

    data = csv.DictReader(
        c.result['Drivers_type_role.tsv'],
        delimiter = '\t',
    )
    result = collections.defaultdict(set)

    for rec in data:

        uniprots = mapping.map_name(
            rec['geneHGNCsymbol'],
            'genesymbol',
            'uniprot',
        )

        for uniprot in uniprots:

            role_prob, curated = (
                (
                    1.0,
                    True,
                )
                if rec['OncodriveROLE_prob'] == 'Manually curated' else
                (
                    common.float_or_nan(rec['OncodriveROLE_prob']),
                    False,
                )
            )

            result[uniprot].add(
                IntogenAnnotation(
                    type = rec['Driver_type'],
                    role = rec['Role'],
                    curated = curated,
                    oncodrive_role_prob = role_prob,
                )
            )

    return dict(result)
