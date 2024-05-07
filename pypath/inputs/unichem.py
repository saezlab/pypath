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
import sys
import textwrap
import collections

import json

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.common as common
import pypath.share.session as session

_logger = session.Logger(name = 'unichem_input')
_log = _logger._log


def unichem_info():
    """
    List of ID types in UniChem. See more details at
    https://www.ebi.ac.uk/unichem/ucquery/listSources.

    Returns
        (list): A list of named tuples, each representing information about
            one UniChem resource.
    """

    UnichemSource = collections.namedtuple(
        'UnichemSource',
        (
            'number',
            'label',
            'name',
            'description',
            'acquisition',
        ),
    )

    url = urls.urls['unichem']['sources']
    c = curl.Curl(url, large = False, silent = False)
    response = json.loads(c.result)

    result = [
        UnichemSource(
            number = s['sourceID'],
            label = s['nameLabel'],
            name = s['name'],
            description = s['description'],
            acquisition = s['lastUpdated'],
        )
        for s in response['sources']
    ]

    return result


def unichem_sources():
    """
    ID type numeric codes and labels in UniChem. For more information see
    `unichem_info`.

    Returns
        (dict): A dict with ID type numeric IDs as keys and ID type labels
            as values.
    """


    return dict(
        (
            s.number,
            s.label,
        )
        for s in unichem_info()
    )


def unichem_mapping(id_type_a, id_type_b):
    """
    Identifier translation data from UniChem.

    Args
        id_type_a (int,str): An ID type in UniChem: either the integer ID or
            the string label of a resource. For a full list see
            `unichem_sources`.
        id_type_b (int,str): An ID type in UniChem, same way as
            `id_type_a`.

    Returns
        (dict): A dictionary with ID translation data, keys are IDs of
            `id_type_a`, values are sets of IDs of `id_type_b`.
    """

    return (
        _unichem_mapping(id_type_a, id_type_b) or
        common.swap_dict(_unichem_mapping(id_type_b, id_type_a))
    )


def _unichem_mapping(id_type_a, id_type_b):
    """
    Identifier translation data from UniChem.

    Args
        id_type_a (int,str): An ID type in UniChem: either the integer ID or
            the string label of a resource. For a full list see
            `unichem_sources`.
        id_type_b (int,str): An ID type in UniChem, same way as
            `id_type_a`.

    Returns
        (dict): A dictionary with ID translation data, keys are IDs of
            `id_type_a`, values are sets of IDs of `id_type_b`.
    """

    src_to_label = unichem_sources()
    label_to_src = common.swap_dict(src_to_label)

    def get_src_id(id_type):

        id_type = str(id_type)
        _id_type = label_to_src.get(id_type, id_type)

        if not _id_type.isdigit() or _id_type not in src_to_label:

            msg = 'No such ID type: `%s`.' % id_type
            _log(msg)
            raise ValueError(msg)

        return _id_type


    id_type_a = get_src_id(id_type_a)
    id_type_b = get_src_id(id_type_b)

    url = urls.urls['unichem']['mapping'] % (id_type_a, id_type_a, id_type_b)
    c = curl.Curl(
        url,
        large = True,
        silent = False,
        slow = True,
        http2 = False,
    )

    if c.status == 404:

        return {}

    result = collections.defaultdict(set)
    _ = next(c.result)

    for r in c.result:

        src_id, tgt_id = r.strip().split('\t')
        result[src_id].add(tgt_id)

    return dict(result)


def info(source):
    """
    Print information about one source.

    Args
        source (int,str): The numeric or string ID of one source.
    """

    source = str(source)

    for s in unichem_info():

        if source in s[:3]:

            info = [
                '%s [number=%s, label=%s]' % (
                    s.name, s.number, s.label,
                ),
                '',
                'Description:',
            ]

            info.extend(textwrap.wrap(s.description))
            info.extend(['', 'Data acquisition:'])
            info.extend(textwrap.wrap(s.acquisition))
            info.append('')

            sys.stdout.write(os.linesep.join(info))
            sys.stdout.flush()

            break
