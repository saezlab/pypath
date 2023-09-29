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
from past.builtins import xrange, range

import os
import re
import gzip
import struct
import collections

try:
    import urllib2
except:
    import urllib.request as urllib2

import pypath.inputs.uniprot_db as uniprot_db
import pypath.share.progress as progress
import pypath.share.curl as curl
import pypath.share.common as common
import pypath.share.cache as cache
import pypath.share.session as session
import pypath.resources.urls as urls
import pypath.utils.taxonomy as taxonomy

_logger = session.Logger(name = 'pfam_input')
_log = _logger._log


def pfam_uniprot(uniprots = None, organism = 9606):
    """
    Mappings between Pfam and UniProt.

    Args
        uniprots (set): The UniProt IDs to query.
        organism (int): NCBI Taxonomy ID of an organism.

    Returns
        A pair of dicts of sets, the first mapping from UniProt ACs to
        Pfam ACs, the second the other way around.
    """

    uniprots = uniprots or uniprot_db.all_swissprots(organism = organism)
    uniprots = common.to_list(uniprots)
    u_pfam = collections.defaultdict(set)
    pfam_u = collections.defaultdict(set)

    if uniprots is not None:

        prg = progress.Progress(
            len(uniprots) / 30,
            'Downloading data from UniProt',
            1,
        )
        data_all = []

        for i in xrange(0, len(uniprots), 200):

            to = i + 200
            thisPart = uniprots[i:to]
            thisPart = ' OR '.join(['accession:%s' % u for u in thisPart])
            get = {
                'query': thisPart,
                'format': 'tab',
                'columns': 'id,database(Pfam)'
            }

            for j in xrange(3):

                c = curl.Curl(urls.urls['uniprot_basic']['url'], get = get)
                data = c.result

                if data:

                    break

            data = data.split('\n')
            del data[0]
            del data[-1]
            data_all += data
            prg.step()

        prg.terminate()

    else:

        organism = taxonomy.ensure_ncbi_tax_id(organism)

        if not organism:

            return None, None

        organismQuery = 'organism:%u AND reviewed:yes' % organism
        get = {
            'query': organismQuery,
            'format': 'tab',
            'columns': 'id,database(Pfam)'
        }

        for j in xrange(3):

            c = curl.Curl(
                urls.urls['uniprot_basic']['url'],
                get = get,
                silent = False,
                outf = 'uniprot-pfam-%u.tab' % organism,
            )
            data_all = c.result

            if data_all:

                break

        data_all = data_all.split('\n')
        del data_all[0]

    for l in data_all:

        l = l.split('\t')

        pfams = re.sub(';$', '', l[1]).strip()
        pfams = common.to_set(pfams.split(';') if pfams else set())
        u_pfam[l[0]].update(pfams)

        for pfam in pfams:

            pfam_u[pfam].add(l[0])

    return dict(u_pfam), dict(pfam_u)


def pfam_regions(
        uniprots = None,
        pfams = None,
        organism = 9606,
        keepfile = True,
        value = 'both',
    ):
    """
    Args
        uniprots (set): UniProt IDs to include in the result. If neither
            this or ``pfams`` provided, all SwissProts for the given
            organism will be queried.
        pfams (set): Pfam IDs to include in the result.
        organism (int): NCBI Taxonomy ID (or any other name) of the organism.
        keepfile (bool):
            Keep the downloaded file in the cache directory.
        value (str): The return value: either "uniprot", "pfam" or "both".
            This is the direction of the mapping "uniprot" returns a dict
            with UniProt IDs as keys, "pfam" the other way around, a dict
            with Pfam IDs as keys, while "both" returns both dicts as a
            tuple.
    """

    url = urls.urls['pfam_up']['url']
    outf = common.suffix(url, '/')
    urlmd5 = common.md5(url)

    cachefile = os.path.join(
        cache.get_cachedir(),
        '%s-%s' % (urlmd5, outf),
    )
    u_pfam = {}
    pfam_u = {}
    uniprots = common.to_set(uniprots)
    pfams = common.to_set(pfams)

    if not uniprots and not pfams:

        organism = taxonomy.ensure_ncbi_tax_id(organism)
        uniprots = uniprot_db.all_swissprots(organism = organism)

    if not os.path.exists(cachefile):

        _log('Downloading `%s` to `%s`.' % (url, cachefile))
        urllib2.urlretrieve(url, cachefile)
        _log('Finished downloading `%s` to `%s`.' % (url, cachefile))

    with open(cachefile, 'rb') as f:

        f.seek(-4, 2)
        gzsize = struct.unpack('<I', f.read())[0]
        prg = progress.Progress(gzsize, 'Processing Pfam domains', 11)

    with gzip.open(cachefile, 'r') as f:

        for l in f:

            prg.step(len(l))
            l = l.strip().split()

            if l[0] in uniprots or l[4] in pfams:

                if value in {'uniprot', 'both'}:

                    if l[0] not in u_pfam:
                        u_pfam[l[0]] = {}
                    if l[4] not in u_pfam[l[0]]:
                        u_pfam[l[0]][l[4]] = []
                    u_pfam[l[0]][l[4]].append({
                        'isoform': int(l[1]),
                        'start': int(l[5]),
                        'end': int(l[6])
                    })

                if value in {'pfam', 'both'}:

                    if l[4] not in pfam_u:
                        pfam_u[l[4]] = {}
                    if l[0] not in pfam_u[l[4]]:
                        pfam_u[l[4]][l[0]] = []
                    pfam_u[l[4]][l[0]].append({
                        'isoform': int(l[1]),
                        'start': int(l[5]),
                        'end': int(l[6])
                    })

    prg.terminate()

    if not keepfile:

        os.remove(cachefile)

    if value == 'uniprot':
        return u_pfam
    elif value == 'pfam':
        return pfam_u
    else:
        return u_pfam, pfam_u


def pfam_names():
    """
    Mappings between Pfam accessions and human readable names.

    Returns
        A pair of dictionaries, the first maps from names to accessions,
        the second from accessions to names.
    """

    c = curl.Curl(urls.urls['pfam_pdb']['url'], silent = False)
    data = c.result
    dname_pfam = collections.defaultdict(set)
    pfam_dname = collections.defaultdict(set)
    data = data.strip().split('\n')
    del data[0]

    for l in data:

        l = l.split('\t')

        if len(l) > 5:

            pfam = common.prefix(l[4], '.')
            name = l[5]

            pfam_dname[pfam].add(name)
            dname_pfam[name].add(pfam)

    return dict(dname_pfam), dict(pfam_dname)


def pfam_pdb():
    """
    Mappings between Pfam and PDB.

    Returns
        A pair of dicts of dicts, the first mapping from PDB IDs to
        Pfam ACs, the second the other way around. Each inner dict contains
        sets of domains as values, each domain defined by the PDB chain ID,
        and its start and end positions.
    """

    PfamDomain = collections.namedtuple(
        'PfamDomain',
        (
            'chain',
            'start',
            'end',
        ),
    )

    c = curl.Curl(urls.urls['pfam_pdb']['url'], silent = False)
    data = c.result

    pdb_pfam = collections.defaultdict(dict)
    pfam_pdb = collections.defaultdict(dict)
    data = data.strip().split('\n')[2:]

    for l in data:

        l = l.split('\t')

        if len(l) > 4:

            pfam = common.prefix(l[4], '.')
            pdb = l[0].lower()
            chain = l[1]
            start = int(common.non_digit.sub('', l[2]))
            end = int(common.non_digit.sub('', l[3]))

            domain = PfamDomain(chain, start, end)
            pdb_pfam[pdb][pfam] = domain
            pfam_pdb[pfam][pdb] = domain

    return dict(pdb_pfam), dict(pfam_pdb)


def _pfam_uniprot(uniprots, infile = None):

    result = {}
    url = urls.urls['pfam_up']['url']
    c = curl.Curl(url, large = True, silent = False)

    prg = progress.Progress(len(uniprots), 'Looking up domains', 1)

    for l in c.result:

        l = l.split('\t')

        if l[0] in uniprots:
            prg.step()

            if l[0] not in result:
                result[l[0]] = {}

            if l[4] not in result[l[0]]:
                result[l[0]][l[4]] = []

            result[l[0]][l[4]].append([l[1], l[5], l[6]])

    prg.terminate()

    return result
