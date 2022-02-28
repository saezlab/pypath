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

import time
import json

import pycurl
import requests
import bioservices.biomodels as biom

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.settings as settings
import pypath.share.session as session

_logger = session.Logger(name = 'biomodels_input')
_log = _logger._log

"""
Colon character in IDs must be escaped. Askerisk too?
"""

def get_single_model(model_id):
    """Get single BioModel using bioservices."""
    bm = biom.BioModels()
    model = bm.get_model(model_id)

    return model

def get_all_models():
    """
    Fetch list of available models from API.

    Returns: 
        dict: A dictionary of models with model identifiers as keys and
        model attributes as values.
    """

    bm = biom.BioModels()
    models = bm.get_all_models()

    return models




def _get_biomodels():

    def init_fun(resp_hdr):

        return ['Cookie: access-token=%s' % resp_hdr['token']]


    t = int(time.time() * 1000) - 3600000

    loginurl = urls.urls['biomodels']['login'] % t

    hdrs = [
        'Host: www.intomics.com',
        'X-Requested-With: XMLHttpRequest',
        settings.get('user_agent'),
        'Accept-Language: en-US,en;q=0.5',
        'DNT: 1',
        'Connection: keep-alive',
        'Referer: https://www.intomics.com/inbio/map/',
        'Accept: */*'
    ]

    c0 = curl.Curl(loginurl, silent = False, large = False,
                   cache = False, req_headers = hdrs)

    hdrs = hdrs[:-2]

    hdrs.extend([
        'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8',
        'Upgrade-Insecure-Requests: 1',
        'Accept-Encoding: gzip'
    ])

    # 'Host: www.intomics.com' -H 'User-Agent: Mozilla/5.0 (X11; Linux x86_64; rv:52.0) Gecko/20100101 Firefox/52.0' -H 'Accept: text/html,application/xhtml+xml,application/xml;q = 0.9,*/*;q = 0.8' -H 'Accept-Language: en-US,en;q = 0.5' --compressed -H 'Cookie: access_token = '"$token" -H 'DNT: 1' -H 'Connection: keep-alive' -H 'Upgrade-Insecure-Requests: 1'

    hdrs.append('Cookie: access-token=%s' % json.loads(c0.result)['token'])

    url = urls.urls['biomodels']['url']

    time.sleep(1)

    c2 = curl.Curl(
        url,
        silent = False,
        large = True,
        req_headers = hdrs,
        cache = False,
        compressed = True,
    )

    return c0, c2


def get_biomodels(verbose = 0):

    t = int(time.time() * 1000) - 3600000

    url   = 'https://www.intomics.com/inbio/map/api/'\
            'get_data?file=InBio_Map_core_2016_09_12.tar.gz'
    login = 'https://www.intomics.com/inbio/api/login_guest?ref=&_= %u' % t

    fp_login = open('biomodels.login.tmp', 'wb')
    fp_biomodels = open('biomodels.tmp.tar.gz', 'wb')

    c0 = pycurl.Curl()
    c0.setopt(pycurl.URL, login)
    c0.setopt(pycurl.WRITEFUNCTION, fp_login.write)

    c0.perform()

    fp_login.close()

    with open('biomodels.login.tmp', 'r') as fp:

        token = json.loads(fp.read())['token']

    _log('Token: %s' % token)

    hdrs = ['Cookie: access-token=%s' % token]

    c1 = pycurl.Curl()
    c1.setopt(pycurl.URL, url)
    c1.setopt(pycurl.WRITEFUNCTION, fp_biomodels.write)
    c1.setopt(pycurl.HTTPHEADER, [h.encode('ascii') for h in hdrs])
    c1.setopt(pycurl.VERBOSE, 1)
    c1.setopt(pycurl.DEBUGFUNCTION, print)

    c1.perform()

    fp_biomodels.close()


def get_biomodels_req():

    t = int(time.time() * 1000) - 3600000

    url   = 'https://www.intomics.com/inbio/map/api/'\
            'get_data?file=InBio_Map_core_2016_09_12.tar.gz'
    login = 'https://www.intomics.com/inbio/api/login_guest?ref=&_=%u' % t

    r0 = requests.get(login)
    token = json.loads(r0.text)['token']
    hdrs = {'Cookie': 'access-token=%s' % token}

    with open('biomodels.tmp.tar.gz', 'wb') as fp:

        r1 = requests.get(url, headers = hdrs, stream = True)

        for block in r1.iter_content(4096):

            fp.write(block)
