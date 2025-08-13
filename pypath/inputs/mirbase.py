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

from typing import Literal

import re
import collections
import itertools

import pypath.inputs.uniprot as uniprot_input
import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.share.common as common
import pypath_common._constants as _const


def _mirbase_table(name: str):

    url = urls.urls['mirbase']['root'] % name
    c = curl.Curl(url, silent = False, large = False)

    for row in c.result.strip('<>/p\n')[:-3].split('<br>'):

        yield row.split('\t')


def mirbase_organisms(
        from_: Literal['ncbi', 'mirbase', 'latin'] = 'ncbi',
        to: Literal['ncbi', 'mirbase', 'latin'] = 'mirbase',
    ) -> dict:

    COLS = {
        'mirbase': 1,
        'latin': 3,
        'ncbi': 4,
    }

    result = {}

    ifrom = COLS[from_]
    ito = COLS[to]

    for row in _mirbase_table('mirna_species'):

        if ifrom == 4 or ito == 4:

            if row[4] == r'\N':

                continue

            row[4] = int(row[4])

        result[row[ifrom]] = row[ito]

    return result


def mirbase_taxid(organism = 9606):

    import pypath.utils.taxonomy as taxonomy

    return taxonomy.ensure_mirbase_name(organism)


def _latin_name(organism = 9606):

    import pypath.utils.taxonomy as taxonomy

    return taxonomy.ensure_latin_name(organism)



def _mirbase_of_organism(table, column, organism = None):

    if organism is not None:

        organism = mirbase_taxid(organism)

    for row in _mirbase_table(table):

        if organism and not row[column].startswith(organism):

            continue

        yield row


def mirbase_mirna(organism = None):

    yield from _mirbase_of_organism('mirna', 2, organism)


def mirbase_mirna_mature(organism = None):

    yield from _mirbase_of_organism('mirna_mature', 1, organism)


def mirbase_mirna_pre_mature(organism = None):

    mirna = {
        row[0]: row
        for row in mirbase_mirna(organism)
    }
    mirna_mature = {
        row[0]: row
        for row in mirbase_mirna_mature(organism)
    }

    for row in _mirbase_table('mirna_pre_mature'):

        if row[0] in mirna:

            yield mirna[row[0]], mirna_mature[row[1]]


def mirbase_mature(organism = 9606):

    for row in mirbase_mirna_mature(organism = organism):

        for i in (1, 2):

            if row[i]:

                yield row[3], row[i]


def mirbase_precursor(organism = 9606):

    for row in mirbase_mirna(organism = organism):

        for i in (2, 3):

            if row[i]:

                yield row[1], row[i]


def mirbase_precursor_to_mature(organism = 9606):

    for row in mirbase_mirna_pre_mature(organism = organism):

        for i, j in itertools.product((2, 3), (1, 2)):

            if row[0][i] and row[1][j]:

                yield row[0][i], row[1][j]


def mirbase_ids(organism = 9606):

    for row in mirbase_mirna_pre_mature(organism = organism):

        yield row[1][3], row[0][1]


def mirbase_mature_all(organism = 9606):

    return [i[3] for i in mirbase_mirna_mature(organism = organism)]


def mirbase_precursor_all(organism = 9606):

    return [i[1] for i in mirbase_mirna(organism = organism)]
