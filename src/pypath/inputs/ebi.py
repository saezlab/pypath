#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2021
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: Dénes Türei (turei.denes@gmail.com)
#           Nicolàs Palacio
#           Olga Ivanova
#           Sebastian Lobentanzer
#           Ahmet Rifaioglu
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

from future.utils import iteritems

import json
import collections

import pypath.share.curl as curl


def ebi_rest(
    url,
    qs = None,
    fields = None,
    paginate = True,
    page_length = 500
):
    """
    Collects data from an EBI REST web service.

    Args:
        url (str): The URL of the web service.
        qs (str,dict): Query string parameters to be appended to the URL.
        fields (list): List of field names to be extracted from the result.
        paginate (bool): Retrieve all pages until the end (if False, only
            one page will be downloaded).
        page_length (int): Number of records on one page.

    Returns:
        A list of named tuples with the requested fields.
    """

    result = []

    qs = qs or {}
    qs['size'] = qs.get('size', page_length)

    record = None

    while True:

        page_url = '%s?%s' % (
            url,
            '&'.join(
                '%s=%s' % it
                for it in qs.items()
            )
        )

        c = curl.Curl(page_url)

        this_result = json.loads(c.result)

        qs['page'] = str(int(this_result['page']['number']) + 1)

        data = next(iteritems(this_result['_embedded']).__iter__())[1]

        for item in data:

            if not record:

                record = collections.namedtuple(
                    'WebserviceRecord',
                    fields or sorted(item.keys())
                )

            result.append(
                record(*(
                    item.get(field, None)
                    for field in record._fields
                ))
            )

        if (
            not paginate or
            int(this_result['page']['number']) + 1 >=
            int(this_result['page']['totalPages'])
        ):

            break

    return result
