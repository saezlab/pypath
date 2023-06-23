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


_logger = session.Logger(name = 'cell_input')
_log = _logger._log


def cell_supplementary(supp_url: str, article_url: str) -> str:
    """
    Downloads a supplementary material from the Cell journal webpage.

    Args
        supp_url:
            URL of the supplementary material.
        article_url:
            URL of the article page.

    Return
        The path of the downloaded file.
    """

    c_nocall = curl.Curl(
        supp_url,
        call = False,
        setup = False,
        process = False,
        silent = True,
    )
    c_nocall.get_cache_file_name()
    path = c_nocall.cache_file_name

    req_headers = []

    if not os.path.exists(path):

        cookies = {}
        init_url = article_url

        for step in range(3):

            c_init = curl.Curl(
                init_url,
                silent = True,
                large = True,
                cache = False,
                follow = False,
                req_headers = req_headers + ['user-agent: curl/7.69.1'],
                bypass_url_encoding = True,
                retries = 1,
                empty_attempt_again = False,
            )

            new_cookies = dict(
                tuple(
                    h.decode().split(':')[1].\
                    split(';')[0].\
                    strip().split('=', maxsplit = 1)
                )
                for h in c_init.resp_headers
                if h.lower().startswith(b'set-cookie')
            )
            cookies.update(new_cookies)
            _ = cookies.pop('__cflb', None)

            for h in c_init.resp_headers:

                if h.lower().startswith(b'location'):

                    init_url = h.decode().split(':', maxsplit = 1)[1].strip()

            req_headers = (
                [
                    'Cookie: %s' % (
                        '; '.join(
                            '%s=%s' % cookie
                            for cookie in iteritems(cookies)
                        )
                    )
                ]
                    if cookies else
                []
            )

            _log(
                'HTTP %u; location: `%s`, cookies: `%s`.' % (
                    c_init.status,
                    init_url,
                    req_headers[0] if req_headers else '',
                )
            )

            if c_init.status != 302:

                break

    c_table = curl.Curl(
        supp_url,
        silent = False,
        large = True,
        empty_attempt_again = False,
        req_headers = req_headers + ['user-agent: curl/7.69.1'],
    )
    path = c_table.cache_file_name
    c_table.fileobj.close()

    return path
