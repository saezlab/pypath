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

import re
import itertools
import collections

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.inputs.common as inputs_common
import pypath.utils.mapping as mapping


def ramilowski_interactions(putative = False):
    """
    Downloads and processes ligand-receptor interactions from
    Supplementary Table 2 of Ramilowski 2015.

    Returns list of lists with ligand and receptor gene symbols, reference
    and resources as elements.
    """

    c = curl.Curl(urls.urls['rami']['url'], silent = False, large = True)
    xlsname = c.fname
    del(c)
    raw = inputs_common.read_xls(xlsname, 'All.Pairs')[1:]

    Ramilowski2015Interaction = collections.namedtuple(
        'Ramilowski2015Interaction',
        (
            'ligand',
            'receptor',
            'references',
            'resources',
        ),
    )

    return [
        Ramilowski2015Interaction(
            ligand = r[1],
            receptor = r[3],
            references = r[13].replace(' ', ''), # references
            resources = ';'.join(
                filter(len, itertools.chain(r[5:11], [r[15]]))
            ),
        )
        for r in raw
        if not r[15].startswith('EXCLUDED') and (
            putative or r[15] != 'putative'
        )
    ]

    return raw


def ramilowski_locations(long_notes = False):
    """
    Subcellular location annotations from Ramilowski 2015.
    """

    reloc = re.compile(
        r'([^\(]+[^\s^\(])'
        r'\s?\('
        r'?(?:(.*[^\)])?)'
        r'\)?'
    )
    resep = re.compile(r'[\.;,]')
    renote = re.compile(r'Note=([- \w\(\),\s\+\./%\'":;]*)')

    sources = (
        (4, 'UniProt'),
        (5, 'HPRD'),
        (7, 'LocTree3'),
        (10, 'Consensus'),
        (11, 'Consensus6'),
    )

    RamilowskiLocation = collections.namedtuple(
        'RamilowskiLocation',
        [
            'location',
            'source',
            'tmh',
            'note',
            'long_note',
        ],
    )

    url = urls.urls['rami']['loc']
    c = curl.Curl(url, silent = False, large = True)

    _ = next(c.result)

    result = collections.defaultdict(set)

    for l in c.result:
        l = l.strip('\n\r').split('\t')

        for idx, source in sources:
            locs = l[idx]

            long_note = None
            mnote = renote.search(locs)

            if mnote:
                long_note = mnote.groups()[0]
                locs = renote.sub('', locs)

            for loc in resep.split(locs):

                if ':' in loc and 'GO:' not in loc:

                    loc = loc.split(':')[-1]

                loc = loc.strip().replace('- ', '-').lower()

                if (
                    not loc or
                    len(loc.split()) > 3 or
                    re.search(r'\d', loc) or
                    loc == 'n/a' or
                    any(
                        w in loc for w in
                        ('tumor',)
                    )
                ):
                    continue

                m = reloc.match(loc)

                if not m:
                    continue

                location, note = m.groups()
                tmh = l[9].strip()

                uniprots = mapping.map_name(l[3], 'uniprot', 'uniprot')

                for uniprot in uniprots:

                    result[uniprot].add(
                        RamilowskiLocation(
                            location = (
                                location.lower().replace('=', '').strip()
                            ),
                            source = source,
                            tmh = int(tmh) if tmh.isdigit() else None,
                            note = note,
                            long_note = long_note if long_notes else None,
                        )
                    )

    return dict(result)
