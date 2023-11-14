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

import json
import math
import collections

from typing import Callable, List, Optional, Union

import glom

import pypath.share.curl as curl
import pypath_common._constants as _const
import pypath.share.session as session
import pypath.inputs.common as inputs_common

_logger = session.Logger(name = 'inputs.ebi')
_log = _logger._log


def ebi_rest(
        url: str,
        qs: Optional[dict] = None,
        fields: Optional[inputs_common.GlomSpec] = None,
        paginate: bool = True,
        page_length: int = 500,
        size_param: str = 'size',
        page_param: str = 'offset',
        by_page: bool = True,
        page_field: inputs_common.GlomSpec = 'page.number',
        total_field: inputs_common.GlomSpec = 'page.totalPages',
        record_name: Optional[str] = None,
    ) -> List[tuple]:
    """
    Collects data from an EBI REST web service.

    Args
        url:
            The URL of the web service.
        qs:
            Query string parameters to be appended to the URL.
        fields:
            Glom spec of the fields to be extracted from the result.
        paginate:
            Retrieve all pages until the end (if False, only
            one page will be downloaded).
        page_length:
            Number of records on one page.
        size_param:
            Query string key for number of records per page.
        page_param:
            Query string key to request a specific page by.
        by_page:
            The pagination works by page numbers (True) or item
            numbers (False).
        page_field:
            Glom spec of the JSON field that contains the current
            page number.
        total_field:
            Glom spec of the JSON field that contains the total
            number of pages.
        record_name:
            Class name for the named tuples in the result.

    Details
        Read more about glom specs here:
        https://glom.readthedocs.io/en/latest/tutorial.html

    Returns
        A list of named tuples with the requested fields.
    """

    result = []

    qs = qs or {}
    qs[size_param] = qs.get(size_param, page_length)
    page = 0
    totalrec = -1

    while True:

        page_url = '%s?%s' % (
            url,
            '&'.join(
                '{}={}'.format(*it)
                for it in qs.items()
            )
        )

        _log(
            'Downloading page %u (total: %s).' % (
                page + 1,
                'unknown'
                    if totalrec < 0 else
                str(math.ceil(totalrec / page_length))
            )
        )

        c = curl.Curl(page_url)
        c.get_headers()
        headers = c.resp_headers_dict

        totalrec = int(headers.get('x-pagination-totalrecords', totalrec))

        if not c.result:

            break

        res = inputs_common.json_read(c.result)

        page = glom.glom(
            res,
            (page_field, int),
            default = page,
        ) + 1

        qs[page_param] = page * (by_page or page_length)

        total = glom.glom(res, (total_field, int), default = 0)

        if fields:

            res = inputs_common.json_extract(c.result, fields)

        res = res if isinstance(res, list) else [res]

        if res == [None] or res == [_const.GLOM_ERROR]:

            break

        result.extend(res)

        if (
            not paginate or
            (total and page > total) or
            (totalrec > 0 and len(result) >= totalrec)
        ):

            break

    record_name = (
        record_name or
        '%sRecord' % url.rsplit('/', maxsplit = 1)[-1].capitalize()
    )

    record = collections.namedtuple(
        record_name,
        sorted({k for i in result for k in i.keys()})
    )

    nested = all(
        isinstance(val, list)
        for section in result
        for val in section.values()
    )

    if not nested:

        result = [
            dict((k, [v]) for k, v in section.items())
            for section in result
        ]

    result = [
        record(*values)
        for section in result
        for values in zip(
            *(section.get(f, None) for f in record._fields)
        )
    ]

    return result
