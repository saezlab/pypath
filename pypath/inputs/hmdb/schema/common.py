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

from typing import Any

import pypath.share.common as common
from pypath.inputs.hmdb import _log


XMLNS = '{http://www.hmdb.ca}'

SIMPLE_FIELDS = {
    'accession',
    'name',
    'version',
}

ARRAY_FIELDS = {
    ('secondary_accessions', 'accession'),
    ('synonyms', 'synonym'),
}

PATHWAYS = (
    'pathways',
    ('pathway', 'findall'),
    {'name', 'smpdb_id', 'kegg_map_id'},
)


def _ref_chunk(container: str = 'references') -> tuple:

    return (
        container,
        ('reference', 'findall'),
        'pubmed_id',
        None,
    )


class Field:

    def __init__(self, name, *definition):

        if isinstance(name, Field):

            name, definition = name.name, name.definition

        self.name = name
        self.d = definition or (name,)


    def process(self, record) -> tuple[tuple[Any]]:

        value = ((record,),)
        name = common.to_tuple(self.name)

        for d in self.d:

            if d == '@':

                value = tuple(
                    v[0]
                        if isinstance(v[0], (list, tuple)) else v
                    for v in value
                )

            elif d == '*' or isinstance(d, tuple):

                if d == '*':

                    d = tuple(sorted(value[0][0].keys()))

                if len(value[0][0]) > 1:

                    _log('List of dicts content encountered.')

                value = value[0][0]
                value = tuple((value.get(k),) for k in d)
                name = tuple(f'{n}__{k}' for n in name for k in d)

            elif isinstance(d, str):

                value = tuple(tuple(v.get(d) for v in vv) for vv in value)

        return name, value


    def __str__(self):

        return self.name
