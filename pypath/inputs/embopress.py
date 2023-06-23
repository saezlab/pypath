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

from typing import List, Optional, Union

import pypath.share.curl as curl
import pypath.inputs.common as inputs_common
import pypath.share.settings as settings
import pypath.share.session as session

_logger = session.Logger(name = 'inputs.embo')
_log = _logger._log


def embopress_supplementary(
        url: str,
        init_url: str,
        sheet: Optional[Union[int, str]] = None,
    ) -> Union[List[List], str]:
    """

    """

    # EMBOpress, why are you so mean??

    req_headers = [settings.get('user_agent')]

    c00 = curl.Curl(url, call = False, process = False)

    if (
        not os.path.exists(c00.cache_file_name) or
        os.path.getsize(c00.cache_file_name) == 0
    ):
        _log('EMBOpress download: requesting website cookie.')

        c0 = curl.Curl(
            init_url,
            silent = True,
            large = False,
            req_headers = req_headers,
            follow = False,
            cache = False,
            alpn = False,
        )

        cookies = []

        if hasattr(c0, 'resp_headers'):

            for hdr in c0.resp_headers:

                if hdr.lower().startswith(b'set-cookie'):

                    cookie = hdr.split(b':')[1].split(b';')[0].strip()

                    if cookie not in cookies:

                        cookies.append(cookie.decode('ascii'))

            cookies = '; '.join(cookies)

            req_headers.append('Cookie: %s' % cookies)

            _log('Response header: %s' % str(c0.resp_headers))
            _log('Cookies: %s' % str(cookies))
            _log('Request header: %s' % str(req_headers))

        os.remove(c00.cache_file_name)

    c = curl.Curl(
        url,
        silent = False,
        large = True,
        req_headers = req_headers,
        alpn = False,
    )

    fname = c.fname
    del c

    return (
        fname
            if sheet is None else
        inputs_common.read_xls(fname, sheet = sheet)
    )
