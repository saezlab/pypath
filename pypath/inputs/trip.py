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
import bs4
import pickle

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.share.cache as cache
import pypath.share.session as session
import pypath.share.common as common

_logger = session.Logger(name = 'trip_input')
_log = _logger._log


def take_a_trip(cachefile = None):
    """
    Downloads TRIP data from webpage and preprocesses it.
    Saves preprocessed data into `cachefile` and next
    time loads from this file.

    :arg cachefile str:
        Path to pickle dump of preprocessed TRIP database. If does not exist
        the database will be downloaded and saved to this file. By default
        the path queried from the ``settings`` module.
    """

    cachefile = cachefile or cache.cache_item('trip_preprocessed')

    if os.path.exists(cachefile):

        _log(
            'Loading preprocessed TRIP database '
            'content from `%s`' % cachefile
        )
        result = pickle.load(open(cachefile, 'rb'))

        return result

    _log('No cache found, downloading and preprocessing TRIP database.')

    result = {'sc': {}, 'cc': {}, 'vvc': {}, 'vtc': {}, 'fc': {}}
    intrs = {}
    titles = {
        'Characterization': 'cc',
        'Screening': 'sc',
        'Validation: In vitro validation': 'vtc',
        'Validation: In vivo validation': 'vvc',
        'Functional consequence': 'fc',
    }

    interactors = {}
    base_url = urls.urls['trip']['base_rescued']
    show_url = urls.urls['trip']['show_rescued']
    c = curl.Curl(base_url)
    mainhtml = c.result
    mainsoup = bs4.BeautifulSoup(mainhtml, 'html.parser')
    trppages = common.flat_list(
        [
            [a.attrs['href'] for a in ul.find_all('a')]
            for ul in mainsoup.
                find('div', id = 'trp_selector').
                find('ul').
                find_all('ul')
        ]
    )

    for trpp in trppages:

        trp = trpp.split('/')[-1]
        trpurl = show_url % trp
        c = curl.Curl(trpurl, silent = False)
        trphtml = c.result
        trpsoup = bs4.BeautifulSoup(trphtml, 'html.parser')
        trp_uniprot = trip_find_uniprot(trpsoup)

        if trp_uniprot is None or len(trp_uniprot) < 6:

            _log('Could not find UniProt for %s' % trp)

        for tab in trpsoup.find_all('th', colspan = ['11', '13']):

            ttl = titles[tab.text.strip()]
            tab = tab.find_parent('table')
            trip_process_table(tab, result[ttl], intrs, trp_uniprot)

    _log('Saving processed TRIP database content to `%s`' % cachefile)
    pickle.dump(result, open(cachefile, 'wb'))

    return result


def trip_process_table(tab, result, intrs, trp_uniprot):
    """
    Processes one HTML table downloaded from TRIP webpage.

    @tab : bs4.element.Tag()
        One table of interactions from TRIP webpage.
    @result : dict
        Dictionary the data should be filled in.
    @intrs : dict
        Dictionary of already converted interactor IDs.
        This serves as a cache so do not need to look up
        the same ID twice.
    @trp_uniprot : str
        UniProt ID of TRP domain containing protein.
    """

    for row in tab.find_all('tr'):

        cells = row.find_all(['td', 'th'])

        if 'th' not in [c.name for c in cells]:

            intr = cells[2].text.strip()

            if intr not in intrs:

                intr_uniprot = trip_get_uniprot(intr)
                intrs[intr] = intr_uniprot

                if intr_uniprot is None or len(intr_uniprot) < 6:
                    _log('Could not find UniProt for %s' % intr)

            else:

                intr_uniprot = intrs[intr]

            if (trp_uniprot, intr_uniprot) not in result:

                result[(trp_uniprot, intr_uniprot)] = []

            result[(trp_uniprot, intr_uniprot)].append(
                [c.text.strip() for c in cells]
            )


def trip_get_uniprot(syn):
    """
    Downloads table from TRIP webpage and UniProt attempts to
    look up the UniProt ID for one synonym.

    @syn : str
        The synonym as shown on TRIP webpage.
    """

    url = urls.urls['trip']['show_rescued'] % syn
    c = curl.Curl(url)

    if c.result:

        soup = bs4.BeautifulSoup(c.result, 'html.parser')

        return trip_find_uniprot(soup)


def trip_find_uniprot(soup):
    """
    Looks up a UniProt name in table downloaded from TRIP
    webpage.

    @soup : bs4.BeautifulSoup
        The `BeautifulSoup` instance returned by
        ``pypath.inputs.trip.trip_get_uniprot``.
    """

    for tr in soup.find_all('div', id = 'tab2')[0].find_all('tr'):

        if (
            tr.find('td') is not None and
            tr.find('td').text.strip() == 'Human'
        ):

            uniprot = tr.find_all('td')[2].text.strip()

            return uniprot

    return None


def trip_process(
        exclude_methods = ['Inference', 'Speculation'],
        predictions = False,
        species = 'Human',
        strict = False,
    ):
    """
    Downloads TRIP data by calling `pypath.dadio.take_a_trip()` and
    further provcesses it.
    Returns dict of dict with TRIP data.

    @exclude_methods : list
        Interaction detection methods to be discarded.
    @predictions : bool
        Whether to include predicted interactions.
    @species : str
        Organism name, e.g. `Human`.
    @strict : bool
        Whether include interactions with species not
        used as a bait or not specified.
    """

    nd = 'Not determined'
    spec = set([]) if strict \
        else set(['Not specified', 'Not used as a bait', ''])
    spec.add(species)
    result = {}
    data = take_a_trip()

    for uniprots in common.unique_list(
            common.flat_list([v.keys() for v in data.values()])):
        to_process = False
        refs = set([])
        mets = set([])
        tiss = set([])
        reg = set([])
        eff = set([])

        if uniprots in data['sc']:
            for sc in data['sc'][uniprots]:
                if sc[4] in spec and sc[6] in spec and \
                    (predictions or sc[9] != 'Prediction') and \
                        sc[3] not in exclude_methods:
                    refs.add(sc[10])
                    mets.add(sc[3])
                    tiss.add(sc[7])

        if uniprots in data['vtc']:
            for vtc in data['vtc'][uniprots]:
                if vtc[4] in spec and vtc[7] in spec and \
                        vtc[3] not in exclude_methods:
                    refs.add(vtc[10])
                    mets.add(vtc[3])

        if uniprots in data['vvc']:
            for vvc in data['vvc'][uniprots]:
                if vvc[6] in spec and vvc[8] in spec and \
                        vvc[3] not in exclude_methods:
                    refs.add(vvc[10])
                    mets.add(vvc[3])

                    if len(vvc[4]) > 0:
                        tiss.add(vvc[4])

                    if len(vvc[5]) > 0:
                        tiss.add(vvc[5])

        if uniprots in data['cc']:
            for cc in data['cc'][uniprots]:
                if cc[4] in spec and cc[6] in spec and \
                        cc[3] not in exclude_methods:
                    refs.add(cc[10])
                    mets.add(cc[3])

                    if (cc[5] != nd and len(cc[5]) > 0) or \
                            (cc[7] != nd and len(cc[7]) > 0):
                        reg.add((cc[5], cc[7]))

        if uniprots in data['fc']:
            for fc in data['fc'][uniprots]:
                mets.add(fc[3])
                refs.add(fc[7])

                if len(fc[5]) > 0:
                    eff.add(fc[5])

                if len(fc[6]) > 0:
                    eff.add(fc[6])

        if len(refs) > 0:
            result[uniprots] = {
                'refs': refs,
                'methods': mets,
                'tissues': tiss,
                'effect': eff,
                'regions': reg
            }

    return result


def trip_interactions(
        exclude_methods = ['Inference', 'Speculation'],
        predictions = False,
        species = 'Human',
        strict = False,
    ):
    """
    Obtains processed TRIP interactions by
    calling ``pypath.inputs.trip.trip_process``
    and returns list of interactions. All arguments are passed to
    ``trip_process``, see their definition there.
    """

    data = trip_process(exclude_methods, predictions, species, strict)

    def trip_effect(eff):
        pos = {
            'Sensitization',
            'Activation',
            'Increase in plasma membrane level',
            'Increase in lysosomal membrane level',
            'New channel creation',
        }
        neg = {
            'Desensitization',
            'Decrease in plasma membrane level',
            'Inhibition',
            'Internalization from membrane by ligand',
            'Retain in the endoplasmic reticulum',
        }

        return (
            'stimulation'
                if len(eff & pos) > 0 else
            'inhibition'
                if len(eff & neg) > 0 else
            'unknown'
        )

    return [
        [
            unipr[0],
            unipr[1],
            ';'.join(d['refs']),
            ';'.join(d['methods']),
            trip_effect(d['effect'])
        ]
        for unipr, d in iteritems(data)
    ]
