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

import os

import pypath.share.curl as curl
import pypath.share.session as session


_logger = session.Logger(name = 'science_input')
_log = _logger._log


def science_download(url):
    """
    Downloads a supplementary material from the Science journal webpage.

    Args
        url (str):
            URL of the supplementary material.

    Returns
        The path of the downloaded file.
    """

    c_nocall = curl.Curl(
        url,
        call = False,
        setup = False,
        process = False,
        silent = True,
    )
    c_nocall.get_cache_file_name()
    path = c_nocall.cache_file_name

    req_headers = ['user-agent: curl/7.69.1']

    if not os.path.exists(path):

        c_init = curl.Curl(
            url,
            silent = True,
            large = False,
            cache = False,
            follow = False,
            retries = 1,
            empty_attempt_again = False,
            req_headers = req_headers,
            write_cache = False,
        )

        cookies = dict(
            tuple(
                h.decode().split(':')[1].\
                split(';')[0].\
                strip().split('=', maxsplit = 1)
            )
            for h in c_init.resp_headers
            if h.lower().startswith(b'set-cookie')
        )

        req_headers.append(
            'Cookie: %s' % '; '.join(
                '%s=%s' % ck
                for ck in cookies.items()
            )
        )

        _log(
            'HTTP %u; cookies: `%s`.' % (
                c_init.status,
                req_headers[-1] if req_headers else '',
            )
        )

    c_main = curl.Curl(
        url,
        silent = False,
        large = True,
        empty_attempt_again = False,
        req_headers = req_headers,
    )
    path = c_main.cache_file_name
    c_main.fileobj.close()

    return path
