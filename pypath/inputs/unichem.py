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

import os
import sys
import textwrap
import collections

import bs4

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

    Returns:
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
    soup = bs4.BeautifulSoup(c.result, 'html.parser')
    result = []

    for table in soup.find_all('table'):

        if table.find('tr').text.strip().startswith('src_id'):

            for row in table.find_all('tr')[2:]:

                fields = row.find_all('td')

                result.append(
                    UnichemSource(
                        *(
                            field.text.strip()
                            for field in fields
                        )
                    )
                )

    return result


def unichem_sources():
    """
    ID type numeric codes and labels in UniChem. For more information see
    `unichem_info`.

    Returns:
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


def unichem_mapping(id_type, target_id_type):
    """
    Identifier translation data from UniChem.

    Args:
        id_type (int,str): An ID type in UniChem: either the integer ID or
            the string label of a resource. For a full list see
            `unichem_sources`.
        target_id_type (int,str): An ID type in UniChem, same way as
            `id_type`.

    Returns:
        (dict): A dictionary with ID translation data, keys are IDs of
            `id_type`, values are sets of IDs of `target_id_type`.
    """

    return (
        _unichem_mapping(id_type, target_id_type) or
        common.swap_dict(_unichem_mapping(target_id_type, id_type))
    )


def _unichem_mapping(id_type, target_id_type):
    """
    Identifier translation data from UniChem.

    Args:
        id_type (int,str): An ID type in UniChem: either the integer ID or
            the string label of a resource. For a full list see
            `unichem_sources`.
        target_id_type (int,str): An ID type in UniChem, same way as
            `id_type`.

    Returns:
        (dict): A dictionary with ID translation data, keys are IDs of
            `id_type`, values are sets of IDs of `target_id_type`.
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


    id_type = get_src_id(id_type)
    target_id_type = get_src_id(target_id_type)

    url = urls.urls['unichem']['mapping'] % (id_type, id_type, target_id_type)
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

    Args:
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
