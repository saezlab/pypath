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
import pypath.share.common as common
import pypath.share.session as session


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


class ResidueMapper(session.Logger):
    """
    This class stores and serves the PDB --> UniProt
    residue level mapping. Attempts to download the
    mapping, and stores it for further use. Converts
    PDB residue numbers to the corresponding UniProt ones.
    """


    def __init__(self):

        session.Logger.__init__(self, 'pdb_utils')

        self.clean()


    def load_mapping(self, pdbs):
        """
        Loads PDB-UniProt sequence mapping for one or more PDB IDs.

        Args
            pdb (str,list): One or more PDB IDs.
        """

        non_digit = re.compile(r'[^\d.-]+')
        pdbs = common.to_set(pdbs)
        pdbs = {p.lower() for p in pdbs}

        for pdb in pdbs:

            url = urls.urls['pdb_align']['url'] + pdb

            for attempt in range(3):

                try:

                    data = urllib2.urlopen(url)
                    break

                except:

                    self._log(
                        'Downloading PDB alignment for %s: '
                        '%u attempt failed.' % (pdb, attempt + 1)
                    )

                finally:

                    self._log('Failed to obtain alignment for PDB %s.' % pdb)
                    data = None

            mapper = collections.defaultdict(dict)

            if data:

                alignments = json.loads(data.read())

                for uniprot, alignment in (
                    iteritems(alignments[pdb]['UniProt'])
                ):

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

            self.mappers[pdb] = dict(mapper)


    def get_residue(self, pdb, resnum, chain = None):
        """
        For a residue in a PDB structure returns the UniProt ID and
        the position of the residue in the UniProt sequence.

        Args
            pdb (str): A PDB structure ID.
            resnum (int): The position of the residue.
            chain (str): The chain ID, optional.

        Returns
            Tuple of residue number, offset, UniProt ID and chain ID.
            Returns None if the residue can not be found.
        """

        pdb = pdb.lower()

        if pdb not in self.mappers:

            self.load_mapping(pdb)

        if pdb in self.mappers:

            for _chain, data in iteritems(self.mappers[pdb]):

                pdbends = data.keys()

                if (
                    resnum <= max(pdbends) and (
                        not chain or
                        chain == _chain
                    )
                ):

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

    for l in data:

        l = l.decode('utf-8')

        if not l.startswith('//'):

            l = [x.strip() for x in l.split(':')]
            result[l[0]] = l[1]

    return result