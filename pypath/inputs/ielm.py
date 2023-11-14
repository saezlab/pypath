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

from past.builtins import xrange, range

import os
import sys

try:
    import cPickle as pickle
except:
    import pickle

import bs4

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.progress as progress
import pypath.share.cache as cache_mod
import pypath_common._constants as _const


def get_ielm_huge(
    ppi,
    id_type = 'UniProtKB_AC',
    mydomains = 'HMMS',
    maxwait = 180,
    cache = True,
    part_size = 500,
    headers = None
):
    """
    Loads iELM predicted domain-motif interaction data for a set of
    protein-protein interactions. This method breaks the list into
    reasonable sized chunks and performs multiple requests to iELM,
    and also retries in case of failure, with reducing the request
    size. Provides feedback on the console.

    :param str id_type:
        The type of the IDs in the supplied interaction list.
        Default is 'UniProtKB_AC'.
        Please refer to iELM what type of IDs it does understand.
    :param str mydomains:
        The type of the domain detection method.
        Defaults to 'HMMS'.
        Please refer to iELM for alternatives.
    :param int maxwait:
        The limit of the waiting time in seconds.
    :param bool cache:
        Whether to use the cache or download everything again.
    :param int part_size:
        The number of interactions to be queried in one request.
    :param list headers:
        Additional HTTP headers to send to iELM with each request.
    """

    ranges = range(0, len(ppi), part_size)
    result = []
    done = False

    while not done:

        for r in ranges:

            this_ppi = ppi[r:r + part_size]

            sys.stdout.write(
                '\t:: Part %u/%u: querying %u interactions.\n' % (
                    ranges.index(r) + 1,
                    len(ranges),
                    len(this_ppi),
                )
            )
            sys.stdout.flush()

            this_res = get_ielm(
                this_ppi,
                id_type,
                mydomains,
                maxwait,
                cache,
                part = True,
                headers = headers,
            )

            if this_res:

                if type(this_res) is dict:

                    return this_res

                result += this_res

                if r == ranges[-1]:

                    done = True

            else:
                part_size = max(int(part_size * 0.8), 20)
                ranges = range(r, len(ppi[r:]), part_size)
                sys.stdout.write(
                    '\t:: One query failed. Setting part size to %u\n' %
                    part_size)
                sys.stdout.flush()

                break

    return result


def get_ielm(
    ppi,
    id_type = 'UniProtKB_AC',
    mydomains = 'HMMS',
    maxwait = 180,
    cache = True,
    part = False,
    part_size = 500,
    headers = None
):
    """
    Performs one query to iELM. Parameters are the same as at get_ielm_huge().
    """

    url = urls.urls['proteomic_ielm']['url']
    network = ''
    from_pickle = []
    ppi_pickle = []
    ppi_query = []
    result = []
    pcache = os.path.join(cache_mod.get_cachedir(), 'ielm.pickle')

    if not part and os.path.exists(pcache):

        from_pickle = pickle.load(open(pcache, 'rb'))
        ppi_pickle = from_pickle['ppi']
        ppi_query = list(set(ppi) - set(ppi_pickle))
        result = from_pickle['ielm']

        if len(ppi_query) == 0:
            return result

    else:
        ppi_query = ppi

    if len(ppi_query) > part_size and not part:

        this_result = get_ielm_huge(
            ppi_query,
            id_type,
            mydomains,
            maxwait,
            cache,
            part_size,
            headers,
        )

    for pp in ppi_query:

        network += '%s %s\r\n' % (pp[0], pp[1])

    post = {'network': network, 'databases': id_type, 'mydomains': mydomains}
    net_md5 = common.md5(network)
    cachefile = os.path.join(cache_mod.get_cachedir(), net_md5 + '.ielm')

    if os.path.exists(cachefile) and cache:

        with open(cachefile, 'r') as f:

            data = f.read()

        soup = bs4.BeautifulSoup(data, 'html.parser')
        src = 'cache'

    else:

        c = curl.Curl(
            url, post = post, silent = False, cache = False, req_headers = headers)
        data = c.result
        soup = bs4.BeautifulSoup(data, 'html.parser')
        sessid = soup.find('input', {'name': 'session_ID'})['value']
        src = 'iELM'

    if data is None:

        sys.stdout.write(_const.ERASE_LINE + _const.CURSOR_UP_ONE)
        sys.stdout.write(
            '\t:: Initial query failed. No data retrieved from iELM.\n')
        sys.stdout.flush()

        return None

    wait = 0

    while soup.title.text == 'iELM Wait Page' and wait < maxwait:

        sys.stdout.write(_const.ERASE_LINE + _const.CURSOR_UP_ONE)
        sys.stdout.write('\t:: Waiting for result. Wait time: %u sec. '
                         'Max waiting time: %u sec.\n' % (wait, maxwait))
        sys.stdout.flush()

        post = {
            'session_ID': sessid,
            'database': id_type,
            'number': '',
            'domains': mydomains,
        }

        c = curl.Curl(
            'http://i.elm.eu.org/wait_2/',
            post = post,
            cache = False,
            req_headers = headers,
        )

        data = c.result

        if data is not None:
            soup = bs4.BeautifulSoup(data, 'html.parser')

        time.sleep(3)
        wait += 3

    if len(soup.find_all('table')) == 0:

        sys.stdout.write(_const.ERASE_LINE + _const.CURSOR_UP_ONE)
        sys.stdout.write('\t:: No data retrieved from iELM. \n')
        sys.stdout.flush()
        soup.title.string = 'http://i.elm.eu.org/proteomic_results/%s' % sessid

        return None

    if cache:

        with open(cachefile, 'w') as f:

            f.write(data)

    sys.stdout.write(_const.ERASE_LINE + _const.CURSOR_UP_ONE)
    sys.stdout.write(
        '\t:: Data retrieved from %s in %u seconds.\n' % (src, wait)
    )
    sys.stdout.flush()
    tbl = soup.find('table', {'id': 'example1'})
    this_result = []

    if tbl:

        url = urls.urls['elm_depr']['url']
        depr_c = curl.Curl(url)
        depr_list = depr_c.result
        depr_list = depr_list.replace('"', '').split('\n')[5:]
        depr = [tuple(x.split('\t')) for x in depr_list if len(x) > 0]

        try:
            depr = dict(depr + [tuple([x[0].lower(), x[1]]) for x in depr])

        except:
            print('\n\n\n', depr, '\n\n\n\n')

        rows = tbl.find_all('tr')
        prg = progress.Progress(
            len(rows),
            'Processing data (%u rows)' % (len(rows) - 1),
            3,
        )

        for tr in tbl.find_all('tr'):
            thisRow = [td.text.strip() for td in tr.find_all('td')]

            if len(thisRow) > 15 and not thisRow[0].startswith('Motif'):

                # replacing deprecated ELM names:
                if thisRow[2].lower() in depr:
                    thisRow[2] = depr[thisRow[2].lower()]

                if thisRow[2].lower() in depr:
                    thisRow[2] = depr[thisRow[2].lower()]

                this_result.append(thisRow)

            prg.step()

        prg.terminate()

    if not part:

        result = {
            'ppi': list(set(ppi_pickle + ppi_query)),
            'ielm': result + this_result
        }
        pickle.dump(result, open(pcache, 'wb'))

    return this_result
