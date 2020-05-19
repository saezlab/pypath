#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2020
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#                  Nicolàs Palacio
#                  Olga Ivanova
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

import os
import re
import itertools

import pypath.share.curl as curl
import pypath.share.session as session
import pypath.resources.urls as urls
import pypath.inputs.common as inputs_common

_logger = session.Logger(name = 'inputs.kirouac2010')
_log = _logger._log


def kirouac2010_interactions():
    """
    Returns tuples of ligand-receptor genesymbol pairs.
    """

    rename = re.compile(r'[A-Z]{2}[A-Z0-9][-A-Z0-9]*')
    rerange = re.compile(r'([0-9])-([0-9])')
    reslash = re.compile(r'.*?([A-Z0-9]{1,3}/[/A-Z0-9]+)')


    def get_names(s):

        names = set()
        prev = None

        for n in s.split():

            m = rename.findall(n)

            if m:
                prev = m
                m = reslash.match(n)

                if m:
                    for post in m.groups()[0].split('/'):
                        for pre in prev:
                            names.add('%s%s' % (pre, post))

                else:
                    m = rerange.match(n)

                    if m:
                        intv = m.groups()

                        for post in range(int(intv[0]), int(intv[1]) + 1):
                            for pre in prev:
                                names.add('%s%u' % (pre, post))

                    else:
                        names.update(prev)

            prev = None

        return names


    init_url = urls.urls['kirouac2010']['init_url']
    req_headers = [
        (
            'User-Agent: Mozilla/5.0 (X11; Linux x86_64; rv:68.0) '
            'Gecko/20100101 Firefox/68.0'
        ),
    ]
    url = urls.urls['kirouac2010']['url']

    c00 = curl.Curl(url, call = False, process = False)

    if (
        not os.path.exists(c00.cache_file_name) or
        os.path.getsize(c00.cache_file_name) == 0
    ):
        _log('Kirouac 2010 download: requesting website cookie.')

        c0 = curl.Curl(
            init_url,
            silent = True,
            large = False,
            req_headers = req_headers,
            follow = False,
            cache = False,
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
    )
    xlsname = c.fname
    del(c)
    tbl = inputs_common.read_xls(xlsname, sheet = 'S12')

    result = []

    for r in tbl[2:]:
        namesA = get_names(r[0])
        namesB = get_names(r[1])

        result.extend(list(itertools.product(namesA, namesB)))

    return result
