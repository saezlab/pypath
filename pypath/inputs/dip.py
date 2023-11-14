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

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.inputs.mitab as mitab
import pypath.share.cache as _cache


def dip_interactions(
        url = None,
        organism = 9606,
        core_only = True,
        direct_only = True,
        small_scale_only = True,
    ):

    strDipCore = 'dip-quality-status:core'
    strDirect = 'direct interaction'
    strPhysInt = 'physical interaction'
    strSmallS = 'small scale'
    url = (
        url or
        (urls.urls['dip']['url_rescued'] % ('CR' if core_only else ''))
    )
    c = curl.Curl(url, silent = False, large = True)
    f = c.result
    i = []
    lnum = 0

    for l in f:

        if lnum == 0:
            lnum += 1

            continue

        l = l.replace('\n', '').replace('\r', '')
        l = l.split('\t')
        specA = int(l[9].split(':')[1].split('(')[0])
        specA = 0 if l[9] == '-' else int(l[9].split(':')[1].split('(')[0])
        specB = 0 if l[10] == '-' else int(l[10].split(':')[1].split('(')[0])
        intProp = mitab.mitab_field_list(l[11])
        intProc = mitab.mitab_field_list(l[15])
        dipLinkId = l[13]
        expEv = mitab.mitab_field_list(l[6])
        conf = l[14]

        if organism is None or (specA == organism and specB == organism):

            if (
                (
                    not core_only or
                    strDipCore in conf
                ) and (
                    not direct_only or
                    strDirect in intProp or
                    strPhysInt in intProp
                ) and (
                    not small_scale_only or
                    strSmallS in intProc
                )
            ):

                pm = l[8].replace('pubmed:', '').split('|')
                pm = [p for p in pm if not p.startswith('DIP')]
                l = [l[0], l[1]]
                uniprotA = mitab.mitab_field_uniprot(l[0])
                uniprotB = mitab.mitab_field_uniprot(l[1])

                if uniprotA and uniprotB:

                    i.append([
                        uniprotA,
                        uniprotB,
                        ';'.join(pm),
                        ';'.join(intProp),
                        ';'.join(expEv),
                        dipLinkId,
                    ])

        lnum += 1

    f.close()

    return i


def dip_login(user, passwd):
    """
    This does not work for unknown reasons.

    In addition, the binary_data parameter of Curl().__init__() has been
    changed, below updates are necessary.
    """

    bdr = '---------------------------8945224391427558067125853467'
    useragent = 'Mozilla/5.0 (X11; U; Linux i686; en-US; rv:43.0) '\
        'Gecko/20110304 Firefox/43.0'
    loginfname = os.path.join(_cache.get_cachedir(), 'dip.logindata.tmp')
    url = urls.urls['dip']['login']
    req_hdrs = ['User-Agent: %s' % useragent]
    c = curl.Curl(
        url,
        cache = False,
        write_cache = False,
        req_headers = req_hdrs,
        return_headers = True,
        debug = True)
    res = c.result
    hdr = c.resp_headers
    cookie = hdr['set-cookie'].split(';')[0]
    cookie2 = '%s%u' % (cookie[:-1], int(cookie[-1]) + 1)
    othercookie = 'DIPID=11133%3A'
    req_hdrs = [
        'Content-type: multipart/form-data; '
        'boundary = %s' % bdr,
        'Cookie: %s' % cookie2,
        'Referer: %s' % url,
        'User-Agent: %s' % useragent,
        'Connection: keep-alive',
    ]
    post = {'lgn': '1', 'login': user, 'pass': passwd, 'Login': 'Login'}
    login = '--%s\r\n\r\nContent-Disposition: form-data; name = "lgn"\r\n\r\n1'\
        '\r\n--%s\r\n\r\nContent-Disposition: form-data; name = "login"\r\n\r\n'\
        '%s\r\n--%s\r\n\r\nContent-Disposition: form-data; name = "pass"\r\n\r'\
        '\n%s\r\n--%s\r\n\r\nContent-Disposition: form-data; name = "Login"\r\n'\
        '\r\nLogin\r\n%s--\r\n' % (bdr, bdr, user, bdr, passwd, bdr, bdr)
    # login = login.replace('\r', '')

    with codecs.open(loginfname, encoding = 'ISO-8859-1', mode = 'w') as f:
        f.write(login)

    c = curl.Curl(
        url,
        cache = False,
        write_cache = False,
        follow = True,
        req_headers = req_hdrs,
        timeout = 10,
        binary_data = loginfname,
        return_headers = True,
        debug = True)
    res = c.result
    hdr = c.resp_headers

    return res, hdr
