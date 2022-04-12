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

from future.utils import iteritems
from past.builtins import xrange, range

import sys
import json
import webbrowser

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.common as common
import pypath.share.progress as progress


def open_pubmed(pmid):
    """
    Opens PubMed record in web browser.

    @pmid : str or int
        PubMed ID
    """

    pmid = str(pmid)
    url = urls.urls['pubmed']['url'] % pmid
    webbrowser.open(url)


def only_pmids(idList, strict = True):
    """
    Return elements unchanged which comply with the PubMed ID format,
    and attempts to translate the DOIs and PMC IDs using NCBI
    E-utils.
    Returns list containing only PMIDs.

    @idList : list, str
        List of IDs or one single ID.
    @strict : bool
        Whether keep in the list those IDs which are not PMIDs,
        neither DOIs or PMC IDs or NIH manuscript IDs.
    """
    if type(idList) in common.simple_types:
        idList = [idList]

    pmids = {i for i in idList if isinstance(i, int) or i.isdigit()}
    pmcids = [i for i in idList if i.startswith('PMC')]
    dois = [i for i in idList if '/' in i]
    manuscids = [i for i in idList if i.startswith('NIHMS')]

    if not strict:
        pmids = set(pmids) | set(dois) | set(pmcids) | set(manuscids)

    if len(pmcids) > 0:
        pmids = pmids | set(pmids_list(pmcids))

    if len(dois) > 0:
        pmids = pmids | set(pmids_list(dois))

    return list(pmids)


def get_pmid(idList):
    """
    For a list of doi or PMC IDs
    fetches the corresponding PMIDs.
    """

    if type(idList) in common.simple_types:
        idList = [idList]

    url = urls.urls['pubmed-eutils']['conv'] % ','.join(str(i) for i in idList)
    c = curl.Curl(url, silent = True)
    data = c.result

    try:
        js = json.loads(data)

    except:
        js = {}

    return js


def pmids_dict(idList):

    jsn = get_pmid(idList)
    result = {'doi': {}, 'pmc': {}}

    if 'records' in jsn:
        for r in jsn['records']:
            if 'pmid' in r:
                if 'doi' in r:
                    result['doi'][r['pmid']] = r['doi']

                if 'pmcid' in r:
                    result['pmc'][r['pmid']] = r['pmcid']

    return result


def pmids_list(idList):

    jsn = get_pmid(idList)
    result = []

    if 'records' in jsn:
        for r in jsn['records']:
            if 'pmid' in r:
                result.append(r['pmid'])

    return result


def get_pubmeds(pmids):

    pmids = [str(pmid) for pmid in pmids]
    url = urls.urls['pubmed-eutils']['url']
    cache = len(pmids) < 10
    data = {}
    prg = progress.Progress(
        len(pmids) / 100 + 1,
        'Retrieving data from NCBI e-utils',
        1,
        percent = False)

    for offset in xrange(0, len(pmids), 100):
        prg.step()
        post = {
            'id': ','.join(pmids[offset:offset + 100]),
            'retmode': 'json',
            'db': 'pubmed'
        }

        for i in xrange(3):
            try:
                c = curl.Curl(
                    url,
                    silent = False,
                    cache = cache,
                    post = post,
                    override_post = True,
                )
                res = c.result
                data = dict([(k, v)
                             for k, v in iteritems(json.loads(res)['result'])]
                            + [(k, v) for k, v in iteritems(data)])

                break

            except ValueError:
                sys.stdout.write('\t:: Error in JSON, retry %u\n' % i)
                sys.stdout.flush()

    prg.terminate()

    return data