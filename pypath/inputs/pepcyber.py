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
import collections

import bs4
import pandas as pd

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.progress as progress
import pypath.share.cache as cache_mod
import pypath.share.session as session

_logger = session.Logger(name = 'pepcyber_input')
_log = _logger._log


def pepcyber_interactions(dataframe = False, cache = True):
    """
    Downloads phosphoprotein binding protein interactions
    from the PepCyber database (http://www.pepcyber.org/).

    Args
        dataframe (bool): Return a pandas data frame instead of list of
            tuples.
        cache (bool): Read the data from an intermediate cache file,
            if available.
    """

    PepcyberInteraction = collections.namedtuple(
        'PepcyberInteraction',
        (
            'ppdb_class',
            'ppdb_genesymbol',
            'substrate_genesymbol',
            'binding_seq',
            'binding_pos',
            'all_evidences',
            'n_records',
            'category',
            'substrate_residue',
            'ppdb_uniprot',
            'ppdb_refseq',
            'substrate_uniprot',
            'substrate_refseq',
            'evidence',
            'pmid',
        ),
    )

    def get_cells(row):

        cells = row.find_all('td')

        if len(cells) == 10:
            sp = cells[4].find('span')

            if (
                sp is not None and
                'class' in sp.attrs and
                'sequence' in sp.attrs['class']
            ):

                return cells


    cachefile = os.path.join(
        cache_mod.get_cachedir(),
        'pepcyber_details.tsv',
    )

    if cache and os.path.exists(cachefile):

        _log('Reading data from cache file `%s`.' % cachefile)
        tbl = pd.read_csv(cachefile, sep = '\t', dtype = {'pmid': 'string'})
        result = [
            PepcyberInteraction(
                *(f if pd.notna(f) else None for f in row)
            )
            for row in tbl.itertuples(index = False)
        ]


    else:

        _log('Downloading PepCyber data.')

        url = urls.urls['pepcyber']['rescued']
        # this is huge, takes a few minutes!
        c = curl.Curl(
            url,
            silent = False,
            timeout = 600,
            encoding = 'iso-8859-1',
        )
        data = c.result
        soup = bs4.BeautifulSoup(data, 'html.parser')
        rows = soup.find_all('tr')
        result = []

        prg = progress.Progress(
            len(rows),
            'Retrieving and processing PepCyber data',
            7,
        )

        for row in rows:

            prg.step()
            cells = get_cells(row)

            if cells is None:

                continue

            row_txt = [c.text.strip() for c in cells]

            if len(row_txt) > 9 and row_txt[5].isdigit():

                inum = int(row.find('a')['name'])
                row_txt[9] = (
                    None
                        if 'p' not in row_txt[4] else
                    row_txt[4][row_txt[4].index('p') + 1]
                )

                details = pepcyber_details(inum)

                row_txt.extend(
                    details[row_txt[2]]
                        if row_txt[2] in details else
                    [None, None]
                )
                row_txt.extend(
                    details[row_txt[3]]
                        if row_txt[3] in details else
                    [None, None]
                )

                refs = details['_refs'] or [(None,) * 3]

                for ref in refs:

                    this_record = row_txt[1:] + list(ref[1:])
                    this_record[4] = int(this_record[4])
                    this_record[6] = int(this_record[6])
                    result.append(
                        PepcyberInteraction(*this_record)
                    )

        tbl = pd.DataFrame.from_records(
            result,
            columns = PepcyberInteraction._fields,
        )
        _log('Saving data to `%s`.' % cachefile)
        tbl.to_csv(cachefile, sep = '\t', index = False)

    return tbl if dataframe else result


def pepcyber_details(num):
    """
    Retrieves detailed information about an interaction from the PepCyber
    database.

    Returns
        Dict with gene symbols as keys and lists of length 2 as values,
        with UniProt ID and RefSeq protein ID. A special key `_refs`
        holds a list of dictionaries, each with category, evidence type and
        PubMed reference information.
    """

    PepcyberReference = collections.namedtuple(
        'PepcyberReference',
        ('category', 'evidence', 'reference')
    )

    result = {'_refs': []}
    url = urls.urls['pepcyber']['details_rescued'] % num
    c = curl.Curl(url, encoding = 'iso-8859-1')
    data = c.result

    if data:

        soup = bs4.BeautifulSoup(data, 'html.parser')
        gname = None
        prev = ''

        for td in soup.find_all('td'):

            if prev.startswith('Gene name'):

                gname = td.text.strip().split('(')[0]

            if prev.startswith('RefSeq'):

                refseq = td.text.strip()

            if prev.startswith('SwissProt') and gname is not None:

                swprot = td.text.strip()

                if gname and gname[0] != u'\xce':

                    result[gname] = [swprot, refseq]

                gname = None

            prev = td.text.strip()

        if soup.find(text = 'Records:'):

            refs = (
                soup.find(text = 'Records:').
                parent.
                parent.
                parent.
                next_sibling.
                find('table').
                find_all('tr')
            )[1:]

            result['_refs'] = [
                PepcyberReference(
                    *(
                        td.a.a.text if td.find('a') else td.text
                        for td in tr.find_all('td')
                    )
                )
                for tr in refs
            ]

    return result
