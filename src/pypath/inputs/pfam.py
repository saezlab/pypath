#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#  Helps to translate from the mouse data to human data
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

import pypath.inputs.uniprot as uniprot_input
import pypath.share.progress as progress
import pypath.share.curl as curl
import pypath.share.common as common
import pypath.share.cache as cache
import pypath.share.session as session
import pypath.resources.urls as urls
import pypath.utils.taxonomy as taxonomy

_logger = session.Logger(name = 'uniprot_input')
_log = _logger._log


def pfam_uniprot(uniprots = None, organism = 9606):

    if uniprots is None:

        uniprots = uniprot_input.all_uniprots(
            organism = organism,
            swissprot = True,
        )

    u_pfam = {}
    pfam_u = {}

    if uniprots is not None:

        prg = progress.Progress(
            len(uniprots) / 30,
            'Downloading data from UniProt',
            1,
        )
        data_all = []

        for i in xrange(0, len(uniprots), 30):

            to = i + 30
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
                if data is not None:
                    break
            if data is None:
                return None, None
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
            if data_all is not None:
                break

        if data_all is None:
            return None

        data_all = data_all.split('\n')
        del data_all[0]

    for l in data_all:

        l = l.split('\t')

        pfams = re.sub(';$', '', l[1]).strip()
        pfams = pfams.split(';') if pfams else []

        if l[0] not in u_pfam:

            u_pfam[l[0]] = []

        u_pfam[l[0]] += pfams

        for pfam in pfams:

            if pfam not in pfam_u:
                pfam_u[pfam] = []

            pfam_u[pfam].append(l[0])

    return u_pfam, pfam_u


def pfam_regions(
        uniprots = None,
        pfams = None,
        organism = 9606,
        keepfile = True,
        value = 'both',
    ):
    """
    Args:
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

        organism = taxnomy.ensure_ncbi_tax_id(organism)
        uniprots = uniprot_input.all_swissprots(organism = organism)

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

    Returns:
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

    c = curl.Curl(urls.urls['pfam_pdb']['url'], silent = False)
    data = c.result

    if data is None:

        return None, None

    pdb_pfam = {}
    pfam_pdb = {}
    data = data.replace('\r', '').split('\n')
    del data[0]

    for l in data:

        l = l.split('\t')

        if len(l) > 4:

            pfam = l[4].split('.')[0]
            pdb = l[0].lower()
            chain = l[1]
            start = int(common.non_digit.sub('', l[2]))
            end = int(common.non_digit.sub('', l[3]))
            if pdb not in pdb_pfam:
                pdb_pfam[pdb] = {}
            if pfam not in pfam_pdb:
                pfam_pdb[pfam] = {}
            pdb_pfam[pdb][pfam] = [chain, start, end]
            pfam_pdb[pfam][pdb] = [chain, start, end]

    return pdb_pfam, pfam_pdb


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
