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

import re
import bs4
import warnings

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.session as session
import pypath.share.settings as settings

_logger = session.Logger(name = 'inputs.genecards')
_log = _logger._log


_respace = re.compile(r'\s+')
_summary_sources = {
    'Gene Wiki': 'GeneWiki',
    'UniProtKB/Swiss-Prot': 'UniProt',
}


def genecards_datasheet(gene):
    """
    Retrieves a gene (protein) datasheet from GeneCards.
    Returns HTML as string.

    :param str gene:
        A Gene Symbol or UniProt ID.
    """

    url = urls.urls['genecards']['url'] % gene

    c = curl.Curl(
        url,
        silent = True,
        large = False,
        connect_timeout = settings.get('genecards_datasheet_connect_timeout'),
        # timeout = settings.get('genecards_datasheet_timeout'),
        req_headers = [
            'User-Agent: Mozilla/5.0 (X11; U; Linux i686; en-US; rv:134.0) '
            'Gecko/20110304 Firefox/134.0',
        ],
    )

    if c.status not in {0, 200}:

        _log('Failed to retrieve gene card for ID `%s`.' % gene)

        return None

    return c.result


def genecards_soup(gene):
    """
    Retrieves a gene (protein) datasheet from GeneCards.
    Returns ``bs4.BeautifulSoup`` object.

    :param str gene:
        A Gene Symbol or UniProt ID.
    """

    html = genecards_datasheet(gene)

    if html:

        with warnings.catch_warnings():

            warnings.simplefilter('ignore')
            soup = bs4.BeautifulSoup(html)

        return soup


def genecards_summaries(gene):
    """
    Retrieves the summaries from a GeneCards datasheet. Returns a dict with
    the resource names as keys and the summary texts as values.

    :param str gene:
        A Gene Symbol or UniProt ID.
    """

    result = {}

    soup = genecards_soup(gene)

    if not soup:

        return result

    summaries = soup.select_one('section#summaries')

    if summaries:

        for summary in summaries.select('div.gc-subsection'):

            title = summary.select_one('h3').text.strip('\r\n ')

            if title[:7] in {'No data', 'Additio'}:

                continue

            content = _respace.sub(
                ' ',
                ' '.join(
                    par.text
                    for par in summary.select(':not(:nth-child(1))')
                )
            ).strip('\n\r ')

            for gc_name, name in iteritems(_summary_sources):

                if title.startswith(gc_name):

                    title = name
                    break

                title = title.split(maxsplit = 1)[0]

            if content:

                result[title] = content

    return result
