#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2021
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

from future.utils import iteritems

import re
import json
import collections
import urllib

try:
    import urllib2
except:
    import urllib.request as urllib2

try:
    import urlparse
except:
    import urllib.parse
    urlparse = urllib.parse

import bs4

import pypath.resources.urls as urls


Segment = collections.namedtuple(
    'Segment',
    (
        'uniprot',
        'pdb_start',
        'pdb_end',
        'uniprot_start',
        'uniprot_end',
    ),
)


Residue = collections.namedtuple(
    'Residue',
    (
        'uniprot',
        'chain',
        'resnum',
        'offset',
    ),
)


class ResidueMapper(object):
    """
    This class stores and serves the PDB --> UniProt
    residue level mapping. Attempts to download the
    mapping, and stores it for further use. Converts
    PDB residue numbers to the corresponding UniProt ones.
    """


    def __init__(self):

        self.clean()


    def load_mapping(self, pdb):

        non_digit = re.compile(r'[^\d.-]+')
        pdb = pdb.lower()
        url = urls.urls['pdb_align']['url'] + pdb
        data = urllib2.urlopen(url)
        alignments = json.loads(data.read())
        mapper = collections.defaultdict(dict)

        for uniprot, alignment in iteritems(alignments[pdb]['UniProt']):

            for segment in alignment['mappings']:

                chain = segment['chain_id']
                pdbstart = segment['start']['residue_number']
                pdbend = segment['end']['residue_number']
                uniprotstart = segment['unp_start']
                uniprotend = segment['unp_end']

                if chain not in mapper:

                    mapper[chain] = {}

                mapper[chain][pdbend] = Segment(
                    uniprot = uniprot,
                    pdb_start = pdbstart,
                    pdb_end = pdbend,
                    uniprot_start = uniprotstart,
                    uniprot_end = uniprotend,
                )

        self.mappers[pdb] = mapper


    def get_residue(self, pdb, resnum, chain = None):

        pdb = pdb.lower()

        if pdb not in self.mappers:

            self.load_mapping(pdb)

        if pdb in self.mappers:

            for chain, data in iteritems(self.mappers[pdb]):

                pdbends = data.keys()

                if resnum <= max(pdbends):

                    pdbend = min(
                        [x for x in [e - resnum for e in pdbends]
                         if x >= 0]) + resnum
                    seg = data[pdbend]

                    if seg.pdb_start <= resnum:

                        offset = seg.uniprot_start - seg.pdb_start
                        residue = Residue(
                            resnum = resnum + offset,
                            offset = offset,
                            uniprot = seg.uniprot,
                            chain = chain,
                        )

                        return residue

        return None


    def clean(self):
        """
        Removes cached mappings, freeing up memory.
        """

        self.mappers = {}


def residue_pdb(pdb, chain, residue):

    url = urls.urls['pdbsws']['url']
    params = urlparse.urlencode({
        'plain': 1,
        'qtype': 'pdb',
        'id': pdb,
        'chain': chain,
        'res': residue
    })
    data = urllib2.urlopen(url + "?%s" % params)
    result = {}

    return data

    for l in data:

        if not l.startswith('//'):

            l = [x.strip() for x in l.split(':')]
            result[l[0]] = l[1]

    return result