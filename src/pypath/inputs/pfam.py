#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#  Helps to translate from the mouse data to human data
#
#  Copyright
#  2014-2020
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

import gzip
import re

import pypath.inputs.uniprot as uniprot_input
import pypath.share.progress as progress
import pypath.share.curl as curl
import pypath.share.common as common
import pypath.resources.urls as urls
import pypath.utils.taxonomy as taxonomy


def get_pfam(uniprots = None, organism = 9606):

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


def get_pfam_regions(
        uniprots = [],
        pfams = [],
        keepfile = False,
        dicts = 'both',
    ):

    url = urls.urls['pfam_up']['url']
    outf = url.split('/')[-1]
    urlmd5 = common.md5(url)
    if not os.path.exists(settings.get('cachedir')):
        os.makedirs(settings.get('cachedir'))
    cachefile = os.path.join(settings.get('cachedir'), urlmd5 + '-' + outf)
    u_pfam = {}
    pfam_u = {}
    uniprots = set(uniprots)
    pfams = set(pfams)

    if not os.path.exists(cachefile):
        sys.stdout.write(
            '\t:: Downloading data from %s' %
            url.replace('http://', '').replace('ftp://', '').split('/')[0])
        sys.stdout.flush()
        if hasattr(urllib, 'urlretrieve'):
            urllib.urlretrieve(url, cachefile)
        else:
            urllib.request.urlretrieve(url, cachefile)
        sys.stdout.write('\n')

    with open(cachefile, 'rb') as f:
        f.seek(-4, 2)
        gzsize = struct.unpack('<I', f.read())[0]
        prg = progress.Progress(gzsize, 'Processing Pfam domains', 11)

    with gzip.open(cachefile, 'r') as f: # FIXME: Something went wrong here

        for l in f:

            prg.step(len(l))
            l = l.strip().split()

            if l[0] in uniprots or l[4] in pfams:

                if dicts in ['uniprot', 'both']:

                    if l[0] not in u_pfam:
                        u_pfam[l[0]] = {}
                    if l[4] not in u_pfam[l[0]]:
                        u_pfam[l[0]][l[4]] = []
                    u_pfam[l[0]][l[4]].append({
                        'isoform': int(l[1]),
                        'start': int(l[5]),
                        'end': int(l[6])
                    })

                if dicts in ['pfam', 'both']:

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
    if dicts == 'uniprot':
        return u_pfam
    elif dicts == 'pfam':
        return pfam_u
    else:
        return u_pfam, pfam_u


def get_pfam_names():

    c = curl.Curl(urls.urls['pfam_pdb']['url'], silent = False)
    data = c.result
    if data is None:
        return None, None
    dname_pfam = {}
    pfam_dname = {}
    data = data.replace('\r', '').split('\n')
    del data[0]

    for l in data:

        l = l.split('\t')
        if len(l) > 5:
            pfam = l[4].split('.')[0]
            name = l[5]
            if pfam not in pfam_dname:
                pfam_dname[pfam] = []
            if name not in dname_pfam:
                dname_pfam[name] = []
            pfam_dname[pfam].append(name)
            dname_pfam[name].append(pfam)

    for k, v in iteritems(pfam_dname):
        pfam_dname[k] = list(set(v))
    for k, v in iteritems(dname_pfam):
        dname_pfam[k] = list(set(v))

    return dname_pfam, pfam_dname


def get_pfam_pdb():

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
