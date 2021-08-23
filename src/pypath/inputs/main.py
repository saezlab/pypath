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

#
# This module makes possible
# dynamic data integration, downloads
# files from various resources, in standard
# or non-standard text based and xml formats,
# processes them, can parse html.
#

#TODO make a cleaner separation, remove/refactor repeating code

from __future__ import print_function

from future.utils import iteritems
from past.builtins import xrange, range

import pypath.share.session as session_mod

_logger = session_mod.Logger(name = 'inputs')
_log = _logger._log
_console = _logger._console

import urllib

try:
    import urllib2
except ImportError:
    # this works seemless in Py3:
    import urllib.request
    urllib2 = urllib.request

try:
    import urlparse
except:
    # this works seemless in Py3:
    import urllib.parse
    urlparse = urllib.parse

# Py 2/3
try:
    input = raw_input
except NameError:
    pass

try:
    import cPickle as pickle
except:
    import pickle

try:
    from cStringIO import StringIO
except:
    try:
        from StringIO import StringIO
        from StringIO import StringIO as BytesIO
    except:
        from io import BytesIO
        from io import StringIO

import sys
import os
import io
import re
import time
import itertools
import collections
from collections import Counter
import functools
import gzip
import xlrd
import bs4
from lxml import etree
import time
import copy
import struct
import json
import csv
import pycurl
import webbrowser

try:
    import requests
except ImportError:
    _log('Module `requests` not available.', -1)

import codecs
import base64

try:
    import bioservices
except ImportError:
    _log('Module `bioservices` not available.', -1)

from xlrd import open_workbook
from xlrd.biffh import XLRDError

# from this module

import pypath.utils.mapping as mapping
import pypath.utils.reflists as reflists
import pypath.inputs.uniprot as uniprot_input
import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.share.progress as progress
import pypath.share.common as common
import pypath.internals.intera as intera
import pypath.utils.residues as residues
import pypath.share.settings as settings
import pypath.utils.taxonomy as taxonomy
import pypath.utils.homology as homology_mod
import pypath.inputs.pfam as pfam_input
import pypath.inputs.common as inputs_common
from pypath.resources import data_formats

if 'long' not in __builtins__:
    long = int

if 'unicode' not in __builtins__:
    unicode = str

CURSOR_UP_ONE = '\x1b[1A'
ERASE_LINE = '\x1b[2K'


def get_3dcomplex():
    """
    Downloads and preprocesses data from the 3DComplex database.

    Returns dict of dicts where top level keys are PDB IDs, second level
    keys are pairs of tuples of UniProt IDs and values are list with the
    number of amino acids in contact.
    """

    c = curl.Curl(urls.urls['3dcomplexes_contact']['url'], silent = False)
    contact = c.result
    c = curl.Curl(urls.urls['3dcomplexes_correspondancy']['url'], silent = False)
    corresp = c.result
    u_pdb, pdb_u = get_pdb_chains()

    del u_pdb

    if contact is None or corresp is None or pdb_u is None:
        return None

    contact = contact.split('\n')
    corresp = corresp.split('\n')
    del contact[0]
    corr_dict = {}

    for l in corresp:
        l = l.replace('\r', '').split('\t')

        if len(l) > 2:
            pdb = l[0].split('.')[0]

            if pdb not in corr_dict:
                corr_dict[pdb] = {}

            corr_dict[pdb][l[1]] = l[2]

    compl_dict = {}

    for l in contact:
        l = l.replace('\r', '').split('\t')

        if len(l) > 11 and int(l[11]) == 0 and int(l[10]) == 0:
            compl = l[0]
            pdb = compl.split('_')[0]

            if pdb in corr_dict:
                if l[1] in corr_dict[pdb] and l[2] in corr_dict[pdb]:

                    ch1 = corr_dict[pdb][l[1]]
                    ch2 = corr_dict[pdb][l[2]]

                    if pdb in pdb_u and ch1 in pdb_u[pdb]:
                        up1 = pdb_u[pdb][ch1]['uniprot']

                        if pdb in pdb_u and ch2 in pdb_u[pdb]:
                            up2 = pdb_u[pdb][ch2]['uniprot']

                            if compl not in compl_dict:
                                compl_dict[compl] = {}

                            uniprots = [up1, up2]
                            uniprots.sort()
                            uniprots = tuple(uniprots)

                            if uniprots not in compl_dict[compl]:
                                compl_dict[compl][uniprots] = []

                            compl_dict[compl][uniprots].append(float(l[3]))

    return compl_dict


def _3dcomplex_complexes():
    pass


def get_3dc_ddi():
    c = curl.Curl(urls.urls['3dcomplexes_contact']['url'], silent = False)
    contact = c.result
    c = curl.Curl(urls.urls['3dcomplexes_correspondancy']['url'], silent = False)
    corresp = c.result
    u_pdb, pdb_u = get_pdb_chains()
    del u_pdb

    if contact is None or corresp is None or pdb_u is None:
        return None

    contact = contact.split('\n')
    corresp = corresp.split('\n')
    del contact[0]
    corr_dict = {}
    ddi = []
    uniprots = []

    for l in corresp:
        l = l.replace('\r', '').split('\t')

        if len(l) > 2:
            pdb = l[0].split('.')[0]

            if pdb not in corr_dict:
                corr_dict[pdb] = {}

            corr_dict[pdb][l[1]] = l[2]

    prg = progress.Progress(len(contact), 'Collecting UniProts', 9)

    for l in contact:
        prg.step()
        l = l.replace('\r', '').split('\t')

        if len(l) > 11 and int(l[11]) == 0 and int(l[10]) == 0:
            pdb = l[0].split('_')[0]

            if pdb in corr_dict:
                if l[1] in corr_dict[pdb] and l[2] in corr_dict[pdb]:
                    ch1 = corr_dict[pdb][l[1]]
                    ch2 = corr_dict[pdb][l[2]]

                    if pdb in pdb_u and ch1 in pdb_u[pdb]:
                        up1 = pdb_u[pdb][ch1]['uniprot']

                    if pdb in pdb_u and ch2 in pdb_u[pdb]:
                        up2 = pdb_u[pdb][ch2]['uniprot']

                    uniprots += [up1, up2]

    prg.terminate()
    uniprots = list(set(uniprots))
    u_pfam = pfam_input.get_pfam_regions(uniprots, dicts = 'uniprot')
    prg = progress.Progress(len(contact), 'Processing contact information', 9)

    for l in contact:
        prg.step()
        l = l.replace('\r', '').split('\t')

        if len(l) > 11 and int(l[11]) == 0 and int(l[10]) == 0:
            pdb = l[0].split('_')[0]
            pfams1 = list(set([x.split('.')[0] for x in l[7].split(';')]))
            pfams2 = list(set([x.split('.')[0] for x in l[9].split(';')]))

            if pdb in corr_dict:
                if l[1] in corr_dict[pdb] and l[2] in corr_dict[pdb]:
                    ch1 = corr_dict[pdb][l[1]]
                    ch2 = corr_dict[pdb][l[2]]

                    if pdb in pdb_u and ch1 in pdb_u[pdb]:
                        up1 = pdb_u[pdb][ch1]['uniprot']

                        if pdb in pdb_u and ch2 in pdb_u[pdb]:
                            up2 = pdb_u[pdb][ch2]['uniprot']

                            for pfam1 in pfams1:
                                for pfam2 in pfams2:
                                    pfam1_details = [{
                                        'start': None,
                                        'end': None,
                                        'isoform': 1
                                    }]
                                    pfam2_details = [{
                                        'start': None,
                                        'end': None,
                                        'isoform': 1
                                    }]

                                    if up1 in u_pfam and pfam1 in u_pfam[up1]:
                                        pfam1_details = u_pfam[up1][pfam1]

                                    if up2 in u_pfam and pfam2 in u_pfam[up2]:
                                        pfam2_details = u_pfam[up2][pfam2]

                                    for pfam1_d in pfam1_details:
                                        for pfam2_d in pfam2_details:
                                            dom1 = intera.Domain(
                                                protein = up1,
                                                domain = pfam1,
                                                start = pfam1_d['start'],
                                                end = pfam1_d['end'],
                                                isoform = pfam1_d['isoform'],
                                                chains = {pdb: ch1},
                                            )
                                            dom2 = intera.Domain(
                                                protein = up2,
                                                domain = pfam2,
                                                start = pfam2_d['start'],
                                                end = pfam2_d['end'],
                                                isoform = pfam2_d['isoform'],
                                                chains = {pdb: ch2},
                                            )
                                            dd = intera.DomainDomain(
                                                dom1,
                                                dom2,
                                                pdbs = pdb,
                                                sources = '3DComplex',
                                                contact_residues = float(l[3])
                                            )
                                            ddi.append(dd)

    prg.terminate()

    return ddi


def pisa_bonds(lst, chains):
    non_digit = re.compile(r'[^\d.-]+')
    bonds = []

    for bond in lst.find_all('bond'):
        seqnum1 = int(non_digit.sub('', bond.find('seqnum-1').text))
        seqnum2 = int(non_digit.sub('', bond.find('seqnum-2').text))
        res1 = bond.find('res-1').text
        res1 = res1 if res1 not in common.aaletters else common.aaletters[res1]
        res2 = bond.find('res-2').text
        res2 = res2 if res2 not in common.aaletters else common.aaletters[res2]
        chain1 = bond.find('chain-1').text
        chain2 = bond.find('chain-2').text
        uniprot1 = None if chain1 not in chains else chains[chain1]
        uniprot2 = None if chain2 not in chains else chains[chain2]

        if uniprot1 is not None and uniprot2 is not None:
            bonds.append({
                'chain_1': chain1,
                'uniprot_1': uniprot1,
                'res_1': res1,
                'seqnum_1': seqnum1,
                'chain_2': chain2,
                'uniprot_2': uniprot2,
                'res_2': res2,
                'seqnum_2': seqnum2
            })

    return bonds


def get_pisa(pdblist):
    bond_types = {
        'hbonds': 'h-bonds',
        'sbridges': 'salt-bridges',
        'covbonds': 'cov-bonds',
        'ssbonds': 'ss-bonds'
    }
    interfaces = {}
    cachefile = os.path.join(settings.get('cachedir'), 'pisa.pickle')
    u_pdb, pdb_u = get_pdb_chains()

    if os.path.exists(cachefile):
        try:
            interfaces = pickle.load(open(cachefile, 'rb'))

        except:
            pass

    errors = []
    p = 5
    pdblist = list(set(pdblist) - set(interfaces.keys()))
    prg = progress.Progress(
        len(pdblist) / p,
        'Downloading data from PDBe PISA',
        1,
    )

    for i in xrange(0, len(pdblist), p):
        to = i + p
        thisPart = pdblist[i:to]
        url = urls.urls['pisa_interfaces']['url'] + ','.join(thisPart)
        c = curl.Curl(url, cache = False)
        data = c.result

        if data is None:
            msg = 'Could not download: \n\t\t%s' % url
            errors.append(msg)

            continue

        soup = bs4.BeautifulSoup(data, 'html.parser')
        unmapped_residues = []

        for pdb in soup.find_all('pdb_entry'):
            pdb_id = pdb.find('pdb_code').text.lower()
            interfaces[pdb_id] = {}
            chains = {}
            resconv = ResidueMapper()

            if pdb_id in pdb_u:
                for chain, chain_data in iteritems(pdb_u[pdb_id]):
                    chains[chain] = chain_data['uniprot']

                for interface in pdb.find_all('interface'):
                    for b, t in iteritems(bond_types):
                        lst = interface.find(t)

                        if lst is not None:
                            bonds = pisa_bonds(lst, chains)

                            for bond in bonds:
                                uniprots = (
                                    bond['uniprot_1'],
                                    bond['uniprot_2'],
                                )

                                if uniprots not in interfaces[pdb_id]:
                                    css = common.non_digit.sub(
                                        '', interface.find('css').text)
                                    css = (
                                        None if len(css) == 0 else float(css)
                                    )
                                    area = common.non_digit.sub(
                                        '', interface.find('int_area').text)
                                    area = None if len(area) == 0 else float(
                                        area)
                                    solv_en = common.non_digit.sub(
                                        '',
                                        interface.find('int_solv_en').text
                                    )
                                    solv_en = (
                                        None
                                            if len(solv_en) == 0 else
                                        float(solv_en)
                                    )
                                    stab_en = common.non_digit.sub(
                                        '',
                                        interface.find('stab_en').text
                                    )
                                    stab_en = (
                                        None
                                            if len(stab_en) == 0 else
                                        float(stab_en)
                                    )
                                    interfaces[pdb_id][uniprots] = (
                                        intera.Interface(
                                            uniprots[0],
                                            uniprots[1],
                                            source = 'PISA',
                                            pdb = pdb_id,
                                            css = css,
                                            solv_en = solv_en,
                                            area = area,
                                            stab_en = stab_en,
                                        )
                                    )

                                res1 = resconv.get_residue(
                                    pdb_id,
                                    bond['seqnum_1'],
                                )
                                res2 = resconv.get_residue(
                                    pdb_id,
                                    bond['seqnum_2'],
                                )

                                if (
                                    res1 is not None and
                                    res2 is not None and
                                    res1['uniprot'] == uniprots[0] and
                                    res2['uniprot'] == uniprots[1]
                                ):
                                    interfaces[pdb_id][uniprots].add_residues(
                                        (res1['resnum'], bond['res_1'],
                                         uniprots[0]),
                                        (res2['resnum'], bond['res_2'],
                                         uniprots[1]),
                                        typ = b,
                                    )

                                else:
                                    unmapped_residues.append(
                                        (
                                            pdb_id,
                                            bond['seqnum_1'],
                                            bond['seqnum_2'],
                                            uniprots[0],
                                            uniprots[1],
                                        )
                                    )

        pickle.dump(interfaces, open(cachefile, 'wb'), 2)
        prg.step()

    prg.terminate()

    if len(errors) > 0:
        sys.stdout.write('\t:: Failed to download %u files of total %u:\n\n' %
                         (len(errors), len(lst)))

        for e in errors:
            sys.stdout.write('\t' + e + '\n')

        sys.stdout.flush()

    return interfaces, unmapped_residues


def get_3did_ddi(residues = False, ddi_flat = None, organism = 9606):

    if ddi_flat is None:
        c = curl.Curl(urls.urls['3did_ddi']['url'], silent = False)
        data = c.result
        tmpfile = '3did_flat_tmp'

        if data is None:
            return None

        with open(tmpfile, 'w') as f:
            f.write(data)

        lnum = data.count('\n')
        del data

    else:
        tmpfile = ddi_flat

    u_pfam, pfam_u = pfam_input.get_pfam(organism = organism)
    u_pdb, pdb_u = get_pdb_chains()

    if pfam_u is None or pdb_u is None:
        return None

    ddi = {}
    interfaces = {}
    pdblist = {}
    ddi_collect = False
    con_collect = False
    non_digit = re.compile(r'[^\d.-]+')

    with open(tmpfile, 'r') as f:
        prg = progress.Progress(lnum, 'Reading data', 33)

        for l in f:
            prg.step()

            if l.startswith('#=') and con_collect:
                interfaces[(uniprot1, uniprot2, pdb)].append(this_interface)
                con_collect = False

            if l.startswith('#=ID'):
                # new domain pair: attach previous to results:
                if ddi_collect:
                    for u1 in uniprots1:
                        for u2 in uniprots2:
                            if u1 != u2 and len(pdblist) > 0:
                                if (u1, u2) not in ddi:
                                    ddi[(u1, u2)] = {}

                                if (pfam1, pfam2) not in ddi[(u1, u2)]:
                                    ddi[(u1, u2)][(pfam1, pfam2)] = {
                                        'pdbs': pdblist
                                    }

                    ddi_collect = False

                pdblist = {}
                l = l.split('\t')
                pfam1 = l[3].split('(')[1].split('.')[0]
                pfam2 = l[4].split('.')[0]
                uniprots1 = [] if pfam1 not in pfam_u else pfam_u[pfam1]
                uniprots2 = [] if pfam2 not in pfam_u else pfam_u[pfam2]

                if len(set(uniprots1 + uniprots2)) > 1:
                    ddi_collect = True

            elif l.startswith('#=3D'):
                l = l.split('\t')
                pdb = l[1]
                chain1 = l[2].split(':')[0]
                chain2 = l[3].split(':')[0]

                if pdb in pdb_u and \
                        chain1 in pdb_u[pdb] and \
                        chain2 in pdb_u[pdb]:
                    uniprot1 = pdb_u[pdb][chain1]['uniprot']
                    uniprot2 = pdb_u[pdb][chain2]['uniprot']

                    if uniprot1 != uniprot2:
                        if pdb not in pdblist:
                            pdblist[pdb] = []

                        pdblist[pdb] = common.add_to_list(pdblist[pdb],
                                                        (uniprot1, uniprot2))

                    if residues:
                        #res1 = [int(i) for i in l[2].split(':')[1].split('-')]
                        #res2 = [int(i) for i in l[3].split(':')[1].split('-')]
                        if chain1 != chain2:
                            if pdb_u[pdb][chain1]['offset'] is not None and \
                                    pdb_u[pdb][chain2]['offset'] is not None and \
                                    pdb_u[pdb][chain1]['uniprot'] != \
                                    pdb_u[pdb][chain2]['uniprot']:
                                con_collect = True
                                offset1 = pdb_u[pdb][chain1]['offset']
                                offset2 = pdb_u[pdb][chain2]['offset']
                                this_interface = common.Interface(
                                    uniprot1, uniprot2, source = '3DID', pdb = pdb)

                                if (uniprot1, uniprot2, pdb) not in interfaces:
                                    interfaces[(uniprot1, uniprot2, pdb)] = []

                            else:
                                con_collect = False

            elif not residues or not con_collect:
                continue

            else:
                l = l.split('\t')

                if len(l) > 3:
                    rnum1 = int(non_digit.sub('', l[2])) + offset1
                    rnum2 = int(non_digit.sub('', l[3])) + offset2
                    this_interface.add_residues((rnum1, l[0], uniprot1),
                                                (rnum2, l[1], uniprot2))

        prg.terminate()
        prg = progress.Progress(len(ddi), 'Processing interfaces', 99)

        if residues:
            for u, v1 in iteritems(ddi):
                prg.step()

                for d, v2 in iteritems(v1):
                    for p in v2['pdbs'].keys():
                        if (u[0], u[1], p) in interfaces:
                            ddi[u][d]['interfaces'] = interfaces[(u[0], u[1],
                                                                  p)]

        prg.terminate()

    if ddi_flat is None:
        os.remove(tmpfile)

    if residues:
        return ddi, interfaces

    else:
        return ddi


def get_3did(ddi_flat = None, res = True, organism = 9606, pickl = True):
    resultfile = os.path.join(settings.get('cachedir'), '3did_ddi.pickle')

    if pickl and os.path.exists(resultfile):
        result = pickle.load(open(resultfile, 'rb'))

        if len(result) == 1:
            return result

        else:
            return result[0], result[1]

    if ddi_flat is None:
        c = curl.Curl(urls.urls['3did_ddi']['url'], silent = False)
        data = c.result
        tmpfile = '3did_flat_tmp'

        if data is None:
            return None

        with open(tmpfile, 'w') as f:
            f.write(data)

        lnum = data.count('\n')
        del data

    elif os.path.exists(ddi_flat):
        tmpfile = ddi_flat

    else:
        return None

    u_pdb, pdb_u = get_pdb_chains()
    all_unip = set(uniprot_input.all_uniprots(organism = organism))

    if all_unip is None or pdb_u is None:
        return None

    ddi = []
    interfaces = []
    pdb = pdb_prev = intf = None
    skip = True
    non_digit = re.compile(r'[^\d.-]+')
    rmap = residues.ResidueMapper()

    with open(tmpfile, 'r') as f:
        prg = progress.Progress(
            lnum, 'Processing 3DID domain-domain interactions', 33)

        for l in f:
            prg.step()
            l = l.split('\t')

            if l[0].startswith('#=ID'):
                pfam1 = l[3].split('.')[0][2:]
                pfam2 = l[4].split('.')[0]

            elif l[0].startswith('#=3D'):
                pdb_prev = pdb
                skip = True
                pdb = l[1]
                chain1 = l[2][0]
                chain2 = l[3][0]
                uniprot1 = uniprot2 = None

                if pdb != pdb_prev:
                    rmap.clean()

                if pdb in pdb_u:
                    if chain1 in pdb_u[pdb]:
                        uniprot1 = pdb_u[pdb][chain1]['uniprot']

                    if chain2 in pdb_u[pdb]:
                        uniprot2 = pdb_u[pdb][chain2]['uniprot']

                if uniprot1 is not None and uniprot2 is not None and \
                        uniprot1 in all_unip and uniprot2 in all_unip and \
                        uniprot1 != uniprot2:

                    skip = False

                    if intf is not None:
                        interfaces.append(intf)

                    intf = intera.Interface(uniprot1, uniprot2, '3DID', pdb)
                    u1start = u1end = u2start = u2end = {}

                    if l[2].count('-') == 1:
                        start1 = int(non_digit.sub('', l[2][2:].split('-')[0]))
                        end1 = int(non_digit.sub('', l[2][2:].split('-')[1]))
                        u1start = rmap.pdb2uniprot(pdb, start1, chains = chain1)
                        u1end = rmap.pdb2uniprot(pdb, end1, chains = chain1)

                    if l[3].count('-') == 1:
                        start2 = int(non_digit.sub('', l[3][2:].split('-')[0]))
                        end2 = int(non_digit.sub('', l[3][2:].split('-')[1]))
                        u2start = rmap.pdb2uniprot(pdb, start2, chains = chain2)
                        u2end = rmap.pdb2uniprot(pdb, end2, chains = chain2)

                    u1start = None if len (u1start) == 0 else \
                        u1start[chain1]['resnum']
                    u1end = None if len (u1end) == 0 else \
                        u1end[chain1]['resnum']
                    u2start = None if len (u2start) == 0 else \
                        u2start[chain2]['resnum']
                    u2end = None if len (u2end) == 0 else \
                        u2end[chain2]['resnum']
                    dom1 = intera.Domain(
                        uniprot1,
                        domain = pfam1,
                        start = u1start,
                        end = u1end,
                        isoform = 1)
                    dom2 = intera.Domain(
                        uniprot2,
                        domain = pfam2,
                        start = u2start,
                        end = u2end,
                        isoform = 1)
                    dd = intera.DomainDomain(dom1, dom2, [pdb], '3DID')
                    ddi.append(dd)

            elif not skip and res and not l[0].startswith('//'):
                conv1 = rmap.pdb2uniprot(
                    pdb, int(non_digit.sub('', l[2])), chains = chain1)
                conv2 = rmap.pdb2uniprot(
                    pdb, int(non_digit.sub('', l[3])), chains = chain2)

                if len(conv1) > 0 and len(conv2) > 0:
                    intf.add_residues(
                        (conv1[chain1]['resnum'], l[0], uniprot1),
                        (conv2[chain2]['resnum'], l[1], uniprot2))
        interfaces.append(intf)
        prg.terminate()

    if ddi_flat is None:
        os.remove(tmpfile)

    if res:
        pickle.dump([ddi, interfaces], open(resultfile, 'wb'))

        return ddi, interfaces

    else:
        pickle.dump([ddi], open(resultfile, 'wb'))

        return ddi


def get_3did_dmi(dmi_flat = None):
    resultfile = os.path.join(settings.get('cachedir'), '3did_dmi.pickle')

    if os.path.exists(resultfile):
        return pickle.load(open(resultfile, 'rb'))

    if dmi_flat is None:
        c = curl.Curl(urls.urls['3did_dmi']['url'], silent = False)
        data = c.result
        tmpfile = '3did_dmi_flat_tmp'

        if data is None:
            return None

        with open(tmpfile, 'w') as f:
            f.write(data)

        lnum = data.count('\n')
        del data

    elif os.path.exists(dmi_flat):
        tmpfile = dmi_flat

    else:
        return None

    u_pdb, pdb_u = get_pdb_chains()

    if pdb_u is None:
        return None

    dmi = {}
    non_digit = re.compile(r'[^\d.-]+')
    rmap = residues.ResidueMapper()

    with open(tmpfile, 'r') as f:
        prg = progress.Progress(lnum,
            'Processing 3DID domain-motif interactions', 1)

        for l in f:
            prg.step()
            l = l.strip().split()

            if l[0].startswith('#=ID'):
                domain = l[3]

            if l[0].startswith('#=PT'):
                regex = l[1]

            if l[0].startswith('#=3D'):
                pdb = l[1]
                chain1 = l[2].split(':')[0]
                chain2 = l[3].split(':')[0]

                if l[2].count('-') == 1 and l[3].count('-') == 1:
                    pdb_region1 = [
                        int(non_digit.sub('', x))
                        for x in l[2].split(':')[1].split('-')
                    ]
                    pdb_region2 = [
                        int(non_digit.sub('', x))
                        for x in l[3].split(':')[1].split('-')
                    ]
                    u1start = rmap.pdb2uniprot(
                        pdb, pdb_region1[0], chains = chain1)
                    u1end = rmap.pdb2uniprot(
                        pdb, pdb_region1[1], chains = chain1)
                    u2start = rmap.pdb2uniprot(
                        pdb, pdb_region2[0], chains = chain2)
                    u2end = rmap.pdb2uniprot(
                        pdb, pdb_region2[1], chains = chain2)

                    if len(u1start) != 0 and len(u2start) != 0 and \
                            len(u1end) != 0 and len(u2end) != 0:
                        uniprot_key = (u1start[chain1]['uniprot'],
                                       u2start[chain2]['uniprot'])
                        residue_key = (
                            u1start[chain1]['resnum'], u1end[chain1]['resnum'],
                            u2start[chain2]['resnum'], u2end[chain2]['resnum'])

                        if uniprot_key not in dmi:
                            dmi[uniprot_key] = {}

                        if residue_key not in dmi[uniprot_key]:
                            dmi[uniprot_key][residue_key] = []

                        dmi[uniprot_key][residue_key].append({
                            'pdb': pdb,
                            'regex': regex,
                            'instance': l[4],
                            'domain': domain,
                            'contacts': int(non_digit.sub('', l[5])),
                            'topology': int(non_digit.sub('', l[6]))
                        })

        prg.terminate()

    if dmi_flat is None:
        os.remove(tmpfile)

    if len(rmap.download_errors) > 0:
        sys.stdout.write('Failed to download PDB-UniProt mappings for:\n'
                         '%s\n' % ', '.join(rmap.download_errors))
    pickle.dump(dmi, open(resultfile, 'wb'))

    return dmi


def process_3did_dmi():
    dmi = get_3did_dmi()

    if dmi is None:
        return None

    dname_pfam, pfam_dname = pfam_input.get_pfam_names()
    dname_re = re.compile(r'(.*)(_[A-Z]{3}_)(.*)')
    dmi2 = {}
    prg = progress.Progress(len(dmi), 'Processing data', 11)

    for uniprots, dmis in iteritems(dmi):
        prg.step()

        if uniprots not in dmi2:
            dmi2[uniprots] = []

        for regions, dmi_list in iteritems(dmis):
            new = True

            for dm in dmi_list:
                if new:
                    pfam = None
                    dname = None
                    mname = None
                    name_match = dname_re.match(dm['domain'])

                    if name_match:
                        dname = name_match.groups(0)[0]
                        mname = ''.join(name_match.groups(0)[1:])[1:]

                    if dname in dname_pfam:
                        pfam = dname_pfam[dname][0]

                    domain = pfam if pfam is not None else dname
                    domain_name = 'pfam' if pfam is not None else 'domain_name'
                    dom = intera.Domain(
                        uniprots[0],
                        domain = domain,
                        domain_id_type = domain_name,
                        start = regions[0],
                        end = regions[1])
                    mot = intera.Motif(
                        uniprots[1],
                        regions[2],
                        regions[3],
                        instance = dm['instance'],
                        regex = dm['regex'],
                        motif_name = mname)
                    ptm = intera.Ptm(uniprots[1], motif = mot, source = '3DID')
                    dommot = intera.DomainMotif(dom, ptm, sources = '3DID')
                    new = False

                dommot.add_pdbs(dm['pdb'])

            dmi2[uniprots].append(dommot)

    prg.terminate()

    return dmi2


def get_instruct():
    """
    Instruct contains residue numbers in UniProt sequences, it means
    no further calculations of offsets in chains of PDB structures needed.
    Chains are not given, only a set of PDB structures supporting the
    domain-domain // protein-protein interaction.
    """

    non_digit = re.compile(r'[^\d.-]+')
    c = curl.Curl(urls.urls['instruct_human']['url'], silent = False)
    data = c.result

    if data is None:
        return None

    data = data.replace('\r', '').split('\n')
    del data[0]
    instruct = []

    for l in data:
        l = l.split('\t')

        if len(l) > 12:
            domain1 = l[6]
            domain2 = l[7]
            pdb = l[12].split(';')
            uniprot1 = l[0]
            uniprot2 = l[1]
            seq1 = [[non_digit.sub('', n) for n in s.split(',')]
                    for s in l[10].split(';')]
            seq2 = [[non_digit.sub('', n) for n in s.split(',')]
                    for s in l[11].split(';')]
            instruct.append({
                uniprot1: {
                    'pfam': domain1,
                    'chain': None,
                    'seq': seq1
                },
                uniprot2: {
                    'pfam': domain2,
                    'chain': None,
                    'seq': seq2
                },
                'uniprots': [uniprot1, uniprot2],
                'source': 'Instruct',
                'pdb': pdb,
                'references': l[13].split(';')
            })

    return instruct


def get_instruct_offsets():
    """
    These offsets should be understood as from UniProt to PDB.
    """

    non_digit = re.compile(r'[^\d.-]+')
    c = curl.Curl(urls.urls['instruct_offsets']['url'], silent = False)
    data = c.result

    if data is None:
        return None

    data = data.replace('\r', '').split('\n')
    del data[0]
    offsets = {}

    for l in data:
        l = l.split('\t')

        if len(l) > 2:
            pdb = l[0].lower()
            uniprot = l[1]

            try:
                offset = int(non_digit.sub('', l[2]))
                offsets[(pdb, uniprot)] = offset

            except:
                sys.stdout.write('Error processing line:\n')
                sys.stdout.write(l[2])
                sys.stdout.write('\n')
                sys.stdout.flush()

    return offsets


def get_i3d():
    """
    Interaction3D contains residue numbers in given chains in
    given PDB stuctures, so we need to add an offset to get the residue
    numbers valid for UniProt sequences. Offsets can be obtained from
    Instruct, or from the Pfam PDB-chain-UniProt mapping table.
    """

    dname_pfam, pfam_dname = pfam_input.get_pfam_names()

    if dname_pfam is None:
        sys.stdout.write('\n\t:: Could not get Pfam domain names\n\n')

    non_digit = re.compile(r'[^\d.-]+')
    c = curl.Curl(urls.urls['i3d_human']['url'], silent = False)
    data = c.result

    if data is None:
        return None

    data = data.replace('\r', '').split('\n')
    del data[0]
    i3d = []
    prg = progress.Progress(
        len(data), 'Processing domain-domain interactions', 11)

    for l in data:
        prg.step()
        l = l.split('\t')

        if len(l) > 20:
            domain1 = None if l[13] not in dname_pfam else dname_pfam[l[13]]
            domain2 = None if l[20] not in dname_pfam else dname_pfam[l[20]]
            pdb = l[5]
            uniprot1 = l[0]
            uniprot2 = l[1]
            chain1 = l[7]
            seq1 = [[
                int(non_digit.sub('', l[11])), int(non_digit.sub('', l[12]))
            ]]
            chain2 = l[14]
            seq2 = [[
                int(non_digit.sub('', l[18])), int(non_digit.sub('', l[19]))
            ]]
            i3d.append({
                uniprot1: {
                    'pfam': domain1,
                    'chain': chain1,
                    'seq': seq1
                },
                uniprot2: {
                    'pfam': domain2,
                    'chain': chain2,
                    'seq': seq2
                },
                'uniprots': [uniprot1, uniprot2],
                'source': 'I3D',
                'pdb': [pdb],
                'references': []
            })
    prg.terminate()

    return i3d


def get_switches_elm():
    """
    switches.elm is a resource containing functional switches in molecular regulation,
    in domain-motif level resolution, classified into categories according to their
    mechanism.
    """

    residue = re.compile(r'(^[A-Z])([0-9]+)')
    url = data.formats.urls['switches.elm']['url']
    c = curl.Curl(url, silent = False)
    data = c.result

    if data is None:
        return None

    buff = StringIO()
    buff.write(data)
    cols = {
        'intramol': 3,
        'bindingsite_a': 5,
        'bs_a_start': 6,
        'bs_a_end': 7,
        'uniprot_a': 4,
        'uniprot_b': 8,
        'bindingsite_b': 9,
        'bs_b_start': 10,
        'bs_b_end': 11,
        'affected': 12,
        'type': 13,
        'subtype': 14,
        'mechanism': 15,
        'reversible': 16,
        'outcome': 17,
        'outcomedir': 18,
        'modification': 19,
        'modsites': 20,
        'modifiers': 21,
        'effectors': 22,
        'references': 26
    }
    table = inputs_common.read_table(
        cols = cols,
        fileObject = buff,
        sep2 = subf,
        hdr = 1,
    )
    mod_ont = get_ontology('MOD')

    for l in table:
        if l['modification'].startswith('MOD'):
            if l['modification'] in mod_ont:
                l['modification'] = mod_ont[l['modification']]

        l['references'] = [
            x.replace('PMID:', '').strip() for x in l['references']
        ]
        l['modsites'] = [
            (m.group(2), m.group(1))
            for m in
            [residue.match(s.strip()) for s in l['modsites'].split(';')]
        ]
        l['intramol'] = True if l['intramol'].strip() == 'TRUE' else False
        l['bs_a_start'] = [x.split(';') for x in l['bs_a_start'].strip()]
        l['bs_b_start'] = [x.split(';') for x in l['bs_b_start'].strip()]
        l['bs_a_end'] = [x.split(';') for x in l['bs_a_end'].strip()]
        l['bs_b_end'] = [x.split(';') for x in l['bs_b_end'].strip()]
        l['bindingsite_a'] = [x.split(';') for x in l['bindingsite_a'].strip()]
        l['bindingsite_b'] = [x.split(';') for x in l['bindingsite_b'].strip()]
        l['modifiers'] = [
            x.split(':') for x in l['modifiers'].strip().split(';')
        ]
        bs_a_ids = {}
        bs_b_ids = {}
        mod_ids = {}

        for bs in l['bindingsite_a'].split(';'):
            if ':' in bs:
                bs = bs.split(':')

                if bs[0].lower() not in bs_a_ids:
                    bs_a_ids[bs[0].lower()] = []

                bs_a_ids[bs[0].lower()].append(bs[1])

        for bs in l['bindingsite_b'].split(';'):
            if ':' in bs:
                bs = bs.split(':')

                if bs[0].lower() not in bs_b_ids:
                    bs_b_ids[bs[0].lower()] = []

                bs_b_ids[bs[0].lower()].append(bs[1])

        for mod in l['modifiers'].split(';'):
            if ':' in mod:
                mod = mod.split(':')

                if mod[0].lower() not in mod_ids:
                    mod_ids[mod[0].lower()] = []

                mod_ids[mod[0].lower()].append(mod[1])

        l['bindingsite_a'] = bs_a_ids
        l['bindingsite_b'] = bs_b_ids
        l['modifiers'] = mod_ids

    return table


def get_csa(uniprots = None):
    """
    Downloads and preprocesses catalytic sites data.
    This data tells which residues are involved in the catalytic
    activity of one protein.
    """

    url = urls.urls['catalytic_sites']['url']
    c = curl.Curl(url, silent = False)
    data = c.result

    if data is None:
        return None

    u_pdb, pdb_u = get_pdb_chains()
    buff = StringIO()
    buff.write(data)
    cols = {
        'pdb': 0,
        'id': 1,
        'resname': 2,
        'chain': 3,
        'resnum': 4,
        'chem_fun': 5,
        'evidence': 6,
    }
    table = inputs_common.read_table(
        cols = cols,
        fileObject = buff,
        sep = ',',
        hdr = 1,
    )
    css = {}
    prg = progress.Progress(len(table), 'Processing catalytic sites', 11)

    for l in table:
        if l['pdb'] in pdb_u:
            if l['chain'] in pdb_u[l['pdb']]:
                uniprot = pdb_u[l['pdb']][l['chain']]['uniprot']

                if uniprots is None or uniprot in uniprots:
                    offset = pdb_u[l['pdb']][l['chain']]['offset']

                    if offset is not None:
                        l['resnum'] = int(l['resnum']) + offset

                    else:
                        this_res = residue_pdb(l['pdb'], l['chain'],
                                               l['resnum'])

                        if len(this_res) > 0:
                            l['resnum'] = int(this_res['UPCOUNT'])

                        else:
                            l['resnum'] = None

                    if l['resnum'] is not None:
                        if uniprot not in css:
                            css[uniprot] = {}

                        if l['pdb'] not in css[uniprot]:
                            css[uniprot][l['pdb']] = {}

                        if l['id'] not in css[uniprot][l['pdb']]:
                            css[uniprot][l['pdb']][l['id']] = []

                        css[uniprot][l['pdb']][l['id']].append(
                            intera.Residue(l['resname'], l['resnum'], uniprot))

        prg.step()

    prg.terminate()

    return css


def get_ielm_huge(ppi,
                  id_type = 'UniProtKB_AC',
                  mydomains = 'HMMS',
                  maxwait = 180,
                  cache = True,
                  part_size = 500,
                  headers = None):
    """
    Loads iELM predicted domain-motif interaction data for a set of
    protein-protein interactions. This method breaks the list into
    reasonable sized chunks and performs multiple requests to iELM,
    and also retries in case of failure, with reducing the request
    size. Provides feedback on the console.

    :param str id_type:
        The type of the IDs in the supplied interaction list.
        Default is 'UniProtKB_AC'.
        Please refer to iELM what type of IDs it does understand.
    :param str mydomains:
        The type of the domain detection method.
        Defaults to 'HMMS'.
        Please refer to iELM for alternatives.
    :param int maxwait:
        The limit of the waiting time in seconds.
    :param bool cache:
        Whether to use the cache or download everything again.
    :param int part_size:
        The number of interactions to be queried in one request.
    :param list headers:
        Additional HTTP headers to send to iELM with each request.
    """

    ranges = range(0, len(ppi), part_size)
    result = []
    done = False

    while not done:
        for r in ranges:
            this_ppi = ppi[r:r + part_size]
            sys.stdout.write('\t:: Part %u/%u: querying %u interactions.\n' %
                             (ranges.index(r) + 1, len(ranges), len(this_ppi)))
            sys.stdout.flush()
            this_res = get_ielm(
                this_ppi,
                id_type,
                mydomains,
                maxwait,
                cache,
                part = True,
                headers = headers)

            if this_res:
                if type(this_res) is dict:
                    return this_res

                result += this_res

                if r == ranges[-1]:
                    done = True

            else:
                part_size = max(int(part_size * 0.8), 20)
                ranges = range(r, len(ppi[r:]), part_size)
                sys.stdout.write(
                    '\t:: One query failed. Setting part size to %u\n' %
                    part_size)
                sys.stdout.flush()

                break

    return result


def get_ielm(ppi,
             id_type = 'UniProtKB_AC',
             mydomains = 'HMMS',
             maxwait = 180,
             cache = True,
             part = False,
             part_size = 500,
             headers = None):
    """
    Performs one query to iELM. Parameters are the same as at get_ielm_huge().
    """

    url = urls.urls['proteomic_ielm']['url']
    network = ''
    from_pickle = []
    ppi_pickle = []
    ppi_query = []
    result = []
    pcache = os.path.join(settings.get('cachedir'), 'ielm.pickle')

    if not part and os.path.exists(pcache):
        from_pickle = pickle.load(open(pcache, 'rb'))
        ppi_pickle = from_pickle['ppi']
        ppi_query = list(set(ppi) - set(ppi_pickle))
        result = from_pickle['ielm']

        if len(ppi_query) == 0:
            return result

    else:
        ppi_query = ppi

    if len(ppi_query) > part_size and not part:
        this_result = get_ielm_huge(ppi_query, id_type, mydomains, maxwait,
                                    cache, part_size, headers)

    for pp in ppi_query:
        network += '%s %s\r\n' % (pp[0], pp[1])

    post = {'network': network, 'databases': id_type, 'mydomains': mydomains}
    net_md5 = common.md5(network)
    cachefile = os.path.join(settings.get('cachedir'), net_md5 + '.ielm')

    if os.path.exists(cachefile) and cache:
        with open(cachefile, 'r') as f:
            data = f.read()

        soup = bs4.BeautifulSoup(data, 'html.parser')
        src = 'cache'

    else:
        c = curl.Curl(
            url, post = post, silent = False, cache = False, req_headers = headers)
        data = c.result
        soup = bs4.BeautifulSoup(data, 'html.parser')
        sessid = soup.find('input', {'name': 'session_ID'})['value']
        src = 'iELM'

    if data is None:
        sys.stdout.write(ERASE_LINE + CURSOR_UP_ONE)
        sys.stdout.write(
            '\t:: Initial query failed. No data retrieved from iELM.\n')
        sys.stdout.flush()

        return None

    wait = 0

    while soup.title.text == 'iELM Wait Page' and wait < maxwait:
        # and \
        # len([i for i in soup.find_all('font', {'color': '#FF0000'}) if i.text == \
        #'File no longer available']) == 0:
        sys.stdout.write(ERASE_LINE + CURSOR_UP_ONE)
        sys.stdout.write('\t:: Waiting for result. Wait time: %u sec. '
                         'Max waiting time: %u sec.\n' % (wait, maxwait))
        sys.stdout.flush()
        post = {
            'session_ID': sessid,
            'database': id_type,
            'number': '',
            'domains': mydomains
        }
        c = curl.Curl(
            'http://i.elm.eu.org/wait_2/',
            post = post,
            cache = False,
            req_headers = headers)
        data = c.result

        if data is not None:
            soup = bs4.BeautifulSoup(data, 'html.parser')

        time.sleep(3)
        wait += 3

    if len(soup.find_all('table')) == 0:
        sys.stdout.write(ERASE_LINE + CURSOR_UP_ONE)
        sys.stdout.write('\t:: No data retrieved from iELM. \n')
        sys.stdout.flush()
        soup.title.string = 'http://i.elm.eu.org/proteomic_results/%s' % sessid
        # return {'soup': soup, 'post': urllib.urlencode(post), 'netw':
        # network}

        return None

    if cache:
        with open(cachefile, 'w') as f:
            f.write(data)

    sys.stdout.write(ERASE_LINE + CURSOR_UP_ONE)
    sys.stdout.write('\t:: Data retrieved from %s in %u seconds.\n' %
                     (src, wait))
    sys.stdout.flush()
    tbl = soup.find('table', {'id': 'example1'})
    this_result = []

    if tbl:
        url = urls.urls['elm_depr']['url']
        depr_c = curl.Curl(url)
        depr_list = depr_c.result
        depr_list = depr_list.replace('"', '').split('\n')[5:]
        depr = [tuple(x.split('\t')) for x in depr_list if len(x) > 0]

        try:
            depr = dict(depr + [tuple([x[0].lower(), x[1]]) for x in depr])

        except:
            print('\n\n\n', depr, '\n\n\n\n')

        # redepr = re.compile(r'\b(' + '|'.join(depr.keys()) + r')\b') :(
        rows = tbl.find_all('tr')
        prg = progress.Progress(
            len(rows), 'Processing data (%u rows)' % (len(rows) - 1), 3)

        for tr in tbl.find_all('tr'):
            thisRow = [td.text.strip() for td in tr.find_all('td')]

            if len(thisRow) > 15 and not thisRow[0].startswith('Motif'):
                # replacing deprecated ELM names:
                if thisRow[2].lower() in depr:
                    thisRow[2] = depr[thisRow[2].lower()]

                if thisRow[2].lower() in depr:
                    thisRow[2] = depr[thisRow[2].lower()]
                # thisRow[2] = redepr.sub(lambda x: depr[x.group()],
                # thisRow[2]) :(

                this_result.append(thisRow)

            prg.step()

        prg.terminate()

    if not part:
        result = {
            'ppi': list(set(ppi_pickle + ppi_query)),
            'ielm': result + this_result
        }
        pickle.dump(result, open(pcache, 'wb'))

    return this_result


def pfam_uniprot(uniprots, infile = None):

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


def get_cpdb_ltp():
    return get_cpdb(
        ['HPRD', 'BioGRID', 'PhosphoPOINT', 'MINT', 'BIND', 'IntAct'])


def get_cpdb(exclude = None):
    exclude = set(exclude) if type(exclude) is list else exclude
    result = []
    url = urls.urls['cpdb']['url']
    c = curl.Curl(url, silent = False)
    data = c.result
    data = [
        x.split('\t') for x in data.split('\n')
        if not x.startswith('#') and len(x) > 0
    ]

    for l in data:
        participants = l[2].split(',')

        if len(participants) == 2:
            if not exclude or len(set(l[0].split(',')) - exclude) > 0:
                result.append([participants[0], participants[1], l[0], l[1]])

    return result


def get_dgidb_old():
    """
    Deprecated. Will be removed soon.

    Downloads and processes the list of all human druggable proteins.
    Returns a list of GeneSymbols.
    """

    genesymbols = []
    url = urls.urls['dgidb']['main_url']
    c = curl.Curl(url, silent = False)
    html = c.result
    soup = bs4.BeautifulSoup(html, 'html.parser')
    cats = [
        o.attrs['value']
        for o in soup.find('select', {'id': 'gene_categories'})
        .find_all('option')
    ]

    for cat in cats:
        url = urls.urls['dgidb']['url'] % cat
        c = curl.Curl(url)
        html = c.result
        soup = bs4.BeautifulSoup(html, 'html.parser')
        trs = soup.find('tbody').find_all('tr')
        genesymbols.extend([tr.find('td').text.strip() for tr in trs])

    return mapping.map_names(genesymbols, 'genesymbol', 'uniprot')


def reactome_sbml():
    """
    Downloads Reactome human reactions in SBML format.
    Returns gzip.GzipFile object.
    """

    url = urls.urls['reactome']['sbml']
    c = curl.Curl(url, silent = False, large = True)
    sbml = c.result

    return sbml


def reactome_biopax(organism = 9606, cache = True):
    """
    Downloads Reactome human reactions in SBML format.
    Returns File object.
    """

    organisms = {9606: 'Homo_sapiens'}
    unzipped = os.path.join(
        settings.get('cachedir'),
        'reactome_biopax_%s.owl' % organisms[organism]
    )

    if not os.path.exists(unzipped) or not cache:

        fname = '%s.owl' % organisms[organism]

        url = urls.urls['reactome']['biopax_l3']
        c = curl.Curl(
            url,
            silent = False,
            large = True,
            files_needed = [fname]
        )

        fileobj = c.result[fname]

        with open(unzipped, 'w') as _unzipped:

            while True:

                chunk = fileobj.read(4096)

                if not chunk:
                    break

                _unzipped.write(chunk)

        fileobj.close()

    _unzipped = open(unzipped, 'r')

    return _unzipped


def pid_biopax():
    url = urls.urls['nci-pid']['biopax_l3']
    c = curl.Curl(url, silent = False, large = True)

    return c.result


def panther_biopax():
    url = urls.urls['panther']['biopax_l3']
    c = curl.Curl(url, silent = False, large = True).values()

    return c.result


def acsn_biopax():
    url = urls.urls['acsn']['biopax_l3']
    c = curl.Curl(url, silent = False, large = True)

    return c.result


def reactome_bs():
    sbml = reactome_sbml()
    soup = bs4.BeautifulSoup(sbml.read(), 'html.parser')

    return soup


# Process Reactome BioPAX level 3


def get_soup(elem):

    return bs4.BeautifulSoup(etree.tostring(elem), 'html.parser')


def _bp_collect_resources(elem, tag, restype = None):
    rdfpref = '{http://www.w3.org/1999/02/22-rdf-syntax-ns#}'
    rdfres = '%sresource' % rdfpref

    return [
        x.get(rdfres).replace('#', '') for x in elem.iterfind(tag)
        if rdfres in x.attrib and (restype is None or x.get(rdfres).replace(
            '#', '').startswith(restype))
    ]


def reactions_biopax(biopax_file,
                     organism = 9606,
                     protein_name_type = 'UniProt',
                     clean = True):
    """
    Processes a BioPAX file and extracts binary interactions.
    """

    cachefile = os.path.join(
        settings.get('cachedir'), '%s.processed.pickle' %
            os.path.split(biopax_file.name)[1]
        )

    if os.path.exists(cachefile):
        sys.stdout.write('\t:: Loading already processed data\n')
        sys.stdout.flush()

        return pickle.load(open(cachefile, 'rb'))

    # string constants
    bppref = '{http://www.biopax.org/release/biopax-level3.owl#}'
    rdfpref = '{http://www.w3.org/1999/02/22-rdf-syntax-ns#}'
    rdfid = '%sID' % rdfpref
    rdfab = '%sabout' % rdfpref
    rdfres = '%sresource' % rdfpref
    bpprot = '%sProtein' % bppref
    bpcplx = '%sComplex' % bppref
    bpprre = '%sProteinReference' % bppref
    bpreac = '%sBiochemicalReaction' % bppref
    bpcata = '%sCatalysis' % bppref
    bpctrl = '%sControl' % bppref
    bpcoma = '%sComplexAssembly' % bppref
    bppstp = '%sPathwayStep' % bppref
    bpuxrf = '%sUnificationXref' % bppref
    bpstoi = '%sStoichiometry' % bppref
    bppubr = '%sPublicationXref' % bppref
    bppath = '%sPathway' % bppref
    bpfrfe = '%sFragmentFeature' % bppref
    bpseqi = '%sSequenceInterval' % bppref
    bpseqs = '%sSequenceSite' % bppref
    bpmodf = '%sModificationFeature' % bppref
    bpmodv = '%sSequenceModificationVocabulary' % bppref
    bpmphe = '%smemberPhysicalEntity' % bppref
    bperef = '%sentityReference' % bppref
    bpxref = '%sxref' % bppref
    bpdnam = '%sdisplayName' % bppref
    bprelr = '%sRelationshipXref' % bppref
    bpcsto = '%scomponentStoichiometry' % bppref
    bpstoc = '%sstoichiometricCoefficient' % bppref
    bpphye = '%sphysicalEntity' % bppref
    bpcted = '%scontrolled' % bppref
    bpcter = '%scontroller' % bppref
    bpctyp = '%scontrolType' % bppref
    bpleft = '%sleft' % bppref
    bprgth = '%sright' % bppref
    bpsprc = '%sstepProcess' % bppref
    bpfeat = '%sfeature' % bppref
    bpfelo = '%sfeatureLocation' % bppref
    bpibeg = '%ssequenceIntervalBegin' % bppref
    bpiend = '%ssequenceIntervalEnd' % bppref
    bpseqp = '%ssequencePosition' % bppref
    bpmoty = '%smodificationType' % bppref
    bppcom = '%spathwayComponent' % bppref
    bpterm = '%sterm' % bppref
    bpdb = '%sdb' % bppref
    bpid = '%sid' % bppref
    upStr = 'UniProt'
    modvoc = data_formats.reactome_modifications

    # intermediate results
    proteins = {}
    proteinfamilies = {}
    uniprots = {}
    proteinreferences = {}
    complexes = {}
    complexvariations = {}
    stoichiometries = {}
    reactions = {}
    complexassemblies = {}
    catalyses = {}
    controls = {}
    pathways = {}
    pathwaysteps = {}
    publications = {}
    fragmentfeatures = {}
    sequenceintervals = {}
    sequencesites = {}
    modificationfeatures = {}
    modificationvocabulary = {}
    protein_name_type = protein_name_type.lower()

    # processing the XML
    bpf = reactome_biopax(organism = organism)
    bp_filesize = 0

    if hasattr(biopax_file, 'name') and os.path.exists(biopax_file.name):

        bp_filesize = os.path.getsize(biopax_file.name)

        if biopax_file.mode == 'r':

            biopax_file.close()
            biopax_file = open(biopax_file.name, 'rb')

    elif type(biopax_file) is tarfile.ExFileObject:

        bp_filesize = biopax_file.size

    elif type(biopax_file) is gzip.GzipFile:

        f = open(biopax_file.name, 'rb')
        f.seek(-4, 2)
        bp_filesize = struct.unpack('<I', f.read())[0]
        f.close()

    prg = progress.Progress(bp_filesize, 'Processing BioPAX XML', 1)
    fpos = biopax_file.tell()
    bp = etree.iterparse(biopax_file, events = ('end', ))
    used_elements = []

    try:
        for ev, elem in bp:
            new_fpos = biopax_file.tell()
            prg.step(new_fpos - fpos)
            fpos = new_fpos
            _id = elem.get(rdfid) if rdfid in elem.attrib else elem.get(rdfab)

            # Protein
            if elem.tag == bpprot:
                entref = elem.find(bperef)

                if entref is not None:
                    proteins[_id] = {
                        'protein': entref.get(rdfres).replace('#', ''),
                        'seqfeatures': _bp_collect_resources(elem, bpfeat),
                        'modfeatures': _bp_collect_resources(elem, bpfeat)
                    }

                else:
                    proteinfamilies[_id] = _bp_collect_resources(elem, bpmphe)

            # ProteinReference
            elif elem.tag == bpprre:
                proteinreferences[_id] = _bp_collect_resources(elem, bpxref)

            # UnificationXref
            elif elem.tag == bpuxrf or elem.tag == bprelr:
                db = elem.find(bpdb)

                if db is not None:
                    if elem.find(bpdb).text.lower().startswith(
                            protein_name_type):
                        i = elem.find(bpid)

                        if i is not None:
                            uniprots[_id] = i.text
            # Complex
            elif elem.tag == bpcplx:
                if elem.find(bpcsto) is not None:
                    complexes[_id] = _bp_collect_resources(elem, bpcsto)

                else:
                    complexvariations[_id] = _bp_collect_resources(elem,
                                                                   bpmphe)

            # Stoichiometry
            elif elem.tag == bpstoi:
                stoichiometries[_id] = (elem.find(bpphye).get(rdfres).replace(
                    '#', ''), int(float(elem.find(bpstoc).text)))

            # BiochemicalReaction
            elif elem.tag == bpreac:
                reactions[_id] = {
                    'refs': _bp_collect_resources(elem, bpxref),
                    'left': _bp_collect_resources(elem, bpleft),
                    'right': _bp_collect_resources(elem, bprgth)
                }

            # ComplexAssembly
            elif elem.tag == bpcoma:
                complexassemblies[_id] = {
                    'refs': _bp_collect_resources(elem, bpxref),
                    'left': _bp_collect_resources(elem, bpleft),
                    'right': _bp_collect_resources(elem, bprgth)
                }

            # Catalysis
            elif elem.tag == bpcata:
                cter = elem.find(bpcter)
                cted = elem.find(bpcted)

                if cter is not None and cted is not None:
                    typ = elem.find(bpctyp)
                    catalyses[_id] = {
                        'controller': cter.get(rdfres).replace('#', ''),
                        'controlled': cted.get(rdfres).replace('#', ''),
                        'type': '' if typ is None else typ.text
                    }

            # Control
            elif elem.tag == bpctrl:
                cter = elem.find(bpcter)
                cted = elem.find(bpcted)

                if cter is not None and cted is not None:
                    typ = elem.find(bpctyp)
                    controls[_id] = {
                        'refs': _bp_collect_resources(elem, bpxref),
                        'type': typ.text if typ is not None else '',
                        'controller': cter.get(rdfres).replace('#', ''),
                        'controlled': cted.get(rdfres).replace('#', '')
                    }

            # PathwayStep
            elif elem.tag == bppstp:
                pathwaysteps[_id] = _bp_collect_resources(elem, bppstp)

            # PublicationXref
            elif elem.tag == bppubr:
                pmid = elem.find(bpid)

                if pmid is not None:
                    publications[_id] = pmid.text

            # FragmentFeature
            elif elem.tag == bpfrfe:
                fragmentfeatures[_id] = elem.find(bpfelo).get(rdfres).replace(
                    '#', '')

            # SequenceInterval
            elif elem.tag == bpseqi:
                beg = elem.find(bpibeg)
                end = elem.find(bpiend)
                sequenceintervals[_id] = (
                    beg.get(rdfres).replace('#', '') if beg is not None else
                    None, elem.find(bpiend).get(rdfres).replace('#', '')
                    if end is not None else None)

            # SequenceSite
            elif elem.tag == bpseqs:
                seqp = elem.find(bpseqp)

                if seqp is not None:
                    sequencesites[_id] = int(seqp.text)

            # ModificationFeature
            elif elem.tag == bpmodf:
                felo = elem.find(bpfelo)
                moty = elem.find(bpmoty)

                if felo is not None and moty is not None:
                    modificationfeatures[_id] = (
                        elem.find(bpfelo).get(rdfres).replace('#', ''),
                        elem.find(bpmoty).get(rdfres).replace('#', ''))

            # SequenceModificationVocabulary
            elif elem.tag == bpmodv:
                term = elem.find(bpterm)

                if term is not None:
                    modificationvocabulary[_id] = term.text

            # Pathway
            elif elem.tag == bppath:
                try:
                    pathways[_id] = {
                        'reactions': _bp_collect_resources(elem, bppcom),
                        'pathways': _bp_collect_resources(elem, bppcom)
                    }

                except TypeError:
                    sys.stdout.write('Wrong type at element:\n')
                    sys.stdout.write(etree.tostring(elem))
                    sys.stdout.flush()

            if clean:
                used_elements.append(elem)

                if len(used_elements) > 800:
                    for e in used_elements[:400]:
                        e.clear()

                    used_elements = used_elements[400:]

    except etree.XMLSyntaxError as e:
        sys.stdout.write('\n\tWARNING: XML processing error: %s\n' % str(e))
        sys.stdout.flush()

    prg.terminate()
    del bp
    biopax_file.close()

    # # # # # # # # # # # # # # # # # #
    # from intermediate to final results
    prg = progress.Progress(len(proteins), 'Processing proteins', 11)
    proteins_uniprots = {}

    # return proteinreferences, uniprots
    for pref, protein in iteritems(proteins):
        prg.step()

        if protein['protein'] in proteinreferences:
            for prref in proteinreferences[protein['protein']]:
                if prref in uniprots:
                    proteins_uniprots[pref] = uniprots[prref]

    prg.terminate()
    prg = progress.Progress(len(proteins), 'Processing PTMs', 11)
    proteins_modifications = {}

    for pref, protein in iteritems(proteins):
        prg.step()

        for modf in protein['modfeatures']:
            if modf in modificationfeatures:
                if modificationfeatures[modf][0] in sequencesites:
                    if modificationfeatures[modf][1] in modificationvocabulary:
                        if modificationvocabulary[modificationfeatures[modf][
                                1]] in modvoc:
                            if pref not in proteins_modifications:
                                proteins_modifications[pref] = set([])

                            proteins_modifications[pref].add(
                                (sequencesites[modificationfeatures[modf][0]],
                                 modvoc[modificationvocabulary[
                                     modificationfeatures[modf][1]]][1],
                                 modvoc[modificationvocabulary[
                                     modificationfeatures[modf][1]]][0]))

    prg.terminate()

    # build a uniform dict to handle all protein based entities
    # including complexes and variations/families
    entity_uniprot = {}
    prg = progress.Progress(len(proteins_uniprots), 'Processing proteins', 11)

    for pref, protein in iteritems(proteins_uniprots):
        prg.step()
        entity_uniprot[pref] = [{
            'members': [protein],
            'ptms': {} if protein not in proteins_modifications else {
                protein: proteins_modifications[pref]
            }
        }]

    prg.terminate()
    prg = progress.Progress(
        len(proteinfamilies), 'Processing protein families', 11)

    for pfref, prefs in iteritems(proteinfamilies):
        prg.step()
        entity_uniprot[pfref] = []

        for pref in prefs:
            if pref in proteins_uniprots:
                entity_uniprot[pfref].append({
                    'members': [proteins_uniprots[pref]],
                    'ptms': {} if pref not in proteins_modifications else {
                        proteins_uniprots[pref]: proteins_modifications[pref]
                    }
                })

    prg.terminate()

    # return entity_uniprot, complexes, proteins, proteinreferences, uniprots,
    # proteinfamilies, proteins_uniprots, reactions, controls, catalyses,
    # complexassemblies
    del proteins
    del proteinfamilies
    del proteinreferences

    prg = progress.Progress(len(complexes), 'Processing complexes', 11)

    for cref, cplex in iteritems(complexes):
        prg.step()

        if cref not in entity_uniprot:
            process_complex(0, cref, entity_uniprot, complexes,
                            complexvariations, cplex, stoichiometries)

    prg.terminate()
    del complexes
    del stoichiometries
    del proteins_uniprots

    # return entity_uniprot, proteins, proteinreferences, uniprots, complexes, stoichiometries
    # # #
    prg = progress.Progress(
        len(reactions) + len(complexassemblies), 'Processing reactions', 11)
    reactions_uniprots = \
        process_reactions(reactions, entity_uniprot, publications)
    complexassemblies_uniprots = \
        process_reactions(complexassemblies, entity_uniprot, publications)

    del reactions
    del complexassemblies
    # # #

    prg = progress.Progress(
        len(controls) + len(catalyses), 'Processing controls and catalyses',
        11)
    controls_uniprots = _process_controls(
        dict(itertools.chain.from_iterable(
            iteritems(d) for d in (controls, catalyses)
        )),
        entity_uniprot,
        dict(itertools.chain.from_iterable(
            iteritems(d) for d in (
                reactions_uniprots,
                complexassemblies_uniprots,
            )
        )),
        publications,
    )

    for caref, ca in iteritems(complexassemblies_uniprots):
        controls_uniprots[caref] = {
            'type': 'BINDING',
            'refs':
            [publications[r] for r in ca['refs'] if r in publications],
            'controller': None,
            'controlled': ca
        }

    del entity_uniprot
    pickle.dump(controls_uniprots, open(cachefile, 'wb'))

    # return controls_uniprots, entity_uniprot, proteins, proteinreferences,
    # uniprots, complexes, stoichiometries
    return controls_uniprots


def process_reactions(reactions, entity_uniprot, publications):
    result = {}

    for rref, rea in iteritems(reactions):
        result[rref] = {
            'refs':
            [publications[r] for r in rea['refs'] if r in publications],
            'left':
            [entity_uniprot[l] for l in rea['left'] if l in entity_uniprot],
            'right':
            [entity_uniprot[r] for r in rea['right'] if r in entity_uniprot]
        }

    return result


def _process_controls(controls, entity_uniprot, reactions_uniprots,
                      publications):
    result = {}

    for cref, ctrl in iteritems(controls):
        result[cref] = {
            'type': ctrl['type'],
            'refs':
            [publications[r] for r in ctrl['refs'] if r in publications]
            if 'refs' in ctrl else [],
            'controller': entity_uniprot[ctrl['controller']]
            if ctrl['controller'] in entity_uniprot else None,
            'controlled': reactions_uniprots[ctrl['controlled']]
            if ctrl['controlled'] in reactions_uniprots else None
        }

    return result


def process_complex(depth, cref, entity_uniprot, complexes, complexvariations,
                    cplex, stoichiometries):
    log = open('reactome.log', 'a')
    tabs = '\t' * (depth + 1)
    log.write('%sStarting processing %s, depth = %u\n' %
              (tabs[1:], cref, depth))
    this_cplex = [{'members': [], 'ptms': {}}]
    log.write('%sComplex %s have %u member entities\n' %
              (tabs, cref, len(cplex)))

    for stoi in cplex:
        if stoi in stoichiometries:
            ref, num = stoichiometries[stoi]
            log.write('%sNew member entity: %s, stoichiometric coeff: %u\n' %
                      (tabs, ref, num))

            if ref.startswith('Complex') \
                    and ref not in entity_uniprot:
                if ref in complexes:
                    log.write(
                        '%s%s is a complex with %u subentities, and hasn\'t been processed yet\n'
                        % (tabs, ref, len(complexes[ref])))
                    process_complex(depth + 1, ref, entity_uniprot, complexes,
                                    complexvariations, complexes[ref],
                                    stoichiometries)

                if ref in complexvariations:
                    log.write(
                        '%s%s is a complex group with %u variations, and hasn\'t been processed yet\n'
                        % (tabs, ref, len(complexvariations[ref])))
                    entity_uniprot[ref] = []

                    for mref in complexvariations[ref]:
                        if mref not in entity_uniprot and mref in complexes:
                            log.write(
                                '%s%s is a complex with %u subentities, and hasn\'t been processed yet\n'
                                % (tabs, mref, len(complexes[mref])))
                            process_complex(depth + 1, mref, entity_uniprot,
                                            complexes, complexvariations,
                                            complexes[mref], stoichiometries)

                        if mref in entity_uniprot:
                            log.write(
                                '%s%s is now processed, adding it as an instance of %s\n'
                                % (tabs, mref, ref))
                            entity_uniprot[ref].extend(entity_uniprot[mref])

            if ref in entity_uniprot:
                log.write(
                    '%s%s is an already processed entity, with %u variants and %u members\n'
                    % (tabs, ref, len(entity_uniprot[ref]),
                       len(entity_uniprot[ref][0]['members'])
                       if len(entity_uniprot[ref]) > 0 else 0))
                log.write(
                    '%sNumber of variants after processing %s: %u x %u = %u\n'
                    % (tabs, ref, len(this_cplex), len(entity_uniprot[ref]),
                       len(this_cplex) * len(entity_uniprot[ref])))
                this_cplex_new = []

                for var in this_cplex:
                    i = 0

                    for new_member in entity_uniprot[ref]:
                        var_new = copy.deepcopy(var)
                        var_new['members'].extend(new_member['members'] * num)

                        for u, ptm in iteritems(new_member['ptms']):
                            if u not in var_new['ptms']:
                                var_new['ptms'][u] = set([])
                            var_new['ptms'][u] = var_new['ptms'][
                                u] | new_member['ptms'][u]
                        this_cplex_new.append(var_new)
                        i += 1

                this_cplex = this_cplex_new
                log.write('%sNumber of variants after processing %s: %u\n' %
                          (tabs, ref, len(this_cplex)))
                log.write('%sNumber of members in %s: %u\n' %
                          (tabs, cref, len(this_cplex[0]['members'])
                           if len(this_cplex) > 0 else 0))

            else:
                log.write('%sPermanently missing: %s\n' % (tabs, ref))

    log.write('%sFinished processing %s, found %u variants with %u members\n' %
              (tabs[1:], cref, len(this_cplex), len(this_cplex[0]['members'])
               if len(this_cplex) > 0 else 0))

    if cref not in entity_uniprot:
        entity_uniprot[cref] = []

    entity_uniprot[cref].extend(this_cplex)


def reactome_interactions(cacheFile = None, ask = True, **kwargs):
    """
    Downloads and processes Reactome BioPAX.
    Extracts binary interactions.
    The applied criteria are very stringent, yields very few interactions.
    Requires large free memory, approx. 2G.
    """

    cacheFile = os.path.join(
        settings.get('cachedir'),
        'reactome.interactions.pickle'
    ) if cacheFile is None else cacheFile

    if os.path.exists(cacheFile):
        interactions = pickle.load(open(cacheFile, 'rb'))

    elif ask:

        while True:

            sys.stdout.write(
                '\n\tProcessing Reactome requires huge memory.\n'
                '\tPlease hit `y` if you have at least 2G free memory,\n'
                '\tor `n` to omit Reactome.\n'
                '\tAfter processing once, it will be saved in \n'
                '\t%s, so next time can be loaded quickly.\n\n'
                '\tProcess Reactome now? [y/n]\n' % cacheFile)
            sys.stdout.flush()
            answer = input().lower()

            if answer in {'y', 'n'}:

                break

    else:

        answer = 'y'

    if answer == 'y':
        return get_interactions('reactome', **kwargs)

    else:
        return []


def acsn_interactions_2(**kwargs):
    return get_interactions('acsn', **kwargs)


def pid_interactions(**kwargs):
    return get_interactions('pid', **kwargs)


def panther_interactions(**kwargs):
    return get_interactions('panther', **kwargs)


def get_interactions(source, mandatory_refs = True):
    ctrls = get_controls(source)

    return process_controls(ctrls, mandatory_refs)[0]


def get_controls(source, protein_name_type = None):
    name_types = {
        'acsn': 'HGNC',
        'reactome': 'UniProt',
        'pid': 'UniProt',
        'panther': 'UniProt'
    }

    if protein_name_type is None and source in name_types:
        protein_name_type = name_types[source]
    biopax = globals()['%s_biopax' % source]
    bpfile = biopax()

    if type(bpfile) is list:
        result = {}

        for bpf in bpfile:
            result = dict(
                reactions_biopax(
                    bpf, protein_name_type = protein_name_type).items() +
                result.items())

    else:
        result = reactions_biopax(bpfile, protein_name_type = protein_name_type)

    return result


def process_controls(controls, mandatory_refs = True):
    interactions = set([])
    ptms = []
    regulations = []
    prg = progress.Progress(len(controls), 'Processing interactions', 11)

    for c in controls.values():
        prg.step()

        if len(c['refs']) > 0 or not mandatory_refs:
            if c['controller'] is not None and len(c['controller']) > 0:
                for ctr in c['controller']:
                    if len(common.uniq_list(ctr['members'])) == 1:
                        this_ctr = ctr['members'][0].split('-')[0]
                        ctd = c['controlled']

                        if ctd is not None:
                            # ctd['left'] is not None and ctd['right'] is not
                            # None:

                            for leftInst in itertools.product(*ctd['left']):
                                for rightInst in itertools.product(
                                        *ctd['right']):
                                    lr = common.uniq_list(
                                        common.flat_list([
                                            l['members'] for l in leftInst
                                        ] + [r['members'] for r in rightInst]))

                                    if len(lr) == 1:
                                        this_ctd = lr[0].split('-')[0]
                                        interactions.add((
                                            this_ctr, this_ctd, c['type'],
                                            ';'.join(c['refs'] if len(c[
                                                'refs']) > 0 else ctd['refs']),
                                            'directed'))

                                    else:
                                        modDiff = {}
                                        ptmsLeft = set(
                                            [(ptms[0], ptm)
                                             for l in leftInst
                                             for ptms in l['ptms'].items()
                                             for ptm in ptms[1]])
                                        ptmsRight = set(
                                            [(ptms[0], ptm)
                                             for r in rightInst
                                             for ptms in r['ptms'].items()
                                             for ptm in ptms[1]])
                                        ptmsDiff = ptmsLeft ^ ptmsRight
                                        diffUniProts = common.uniq_list(
                                            [ptm[0] for ptm in ptmsDiff])

                                        if len(diffUniProts) == 1:
                                            this_ctd = diffUniProts[0].split(
                                                '-')[0]
                                            interactions.add(
                                                (this_ctr, this_ctd, c['type'],
                                                 ';'.join(c['refs'] if len(c[
                                                     'refs']) > 0 else ctd[
                                                         'refs']), 'directed'))

                                        else:
                                            lefts = [
                                                set(l['members'])
                                                for l in leftInst
                                            ]
                                            rights = [
                                                set(r['members'])
                                                for r in rightInst
                                            ]
                                            onlyLefts = [
                                                l for l in lefts
                                                if l not in rights
                                            ]
                                            onlyRights = [
                                                r for r in rights
                                                if r not in lefts
                                            ]
                                            diffs = []

                                            for l in onlyLefts:
                                                for r in onlyRights:
                                                    diff = l ^ r
                                                    if len(diff) == 1:
                                                        diffs.append(
                                                            list(diff))
                                            diffs = common.uniq_list(
                                                common.flat_list(diffs))

                                            if len(diffs) == 1:
                                                this_ctd = diffs[0].split('-')[
                                                    0]
                                                interactions.add(
                                                    (this_ctr, this_ctd,
                                                     c['type'],
                                                     ';'.join(c['refs'] if len(
                                                         c['refs']) > 0 else
                                                              ctd['refs']),
                                                     'undirected'))

            # if the controller is unknown
            # and the reaction has only 2 proteins
            # these most probably bind each other
            # to form a complex
            else:
                ctd = c['controlled']
                if ctd is not None:
                    for leftInst in itertools.product(*ctd['left']):
                        for rightInst in itertools.product(*ctd['right']):
                            lr = common.uniq_list(
                                common.flat_list([
                                    l['members'] for l in leftInst
                                ] + [r['members'] for r in rightInst]))

                            if len(lr) == 2:
                                interactions.add(
                                    (lr[0].split('-')[0], lr[1].split('-')[0],
                                     c['type'], ';'.join(ctd['refs'])))

    prg.terminate()

    return list(interactions), ptms, regulations

# Process Reactome SBML

def _reactome_id(obj, attr):
    return _reactome_extract_id(obj.attrs[attr])


def _reactome_extract_id(value):
    return int(value.split('_')[1])


def _reactome_res(obj):
    return _reactome_extract_res(obj.attrs['rdf:resource'])


def _reactome_extract_res(value):
    return value.split(':')[-1]


def _reactome_reactions():
    species = {}
    compartments = {}
    reactions = {}
    soup = reactome_bs()
    m = soup.find('model')

    for cp in m.find('listofcompartments').find_all('compartment'):
        compartments[_reactome_id(cp, 'id')] = cp.attrs['name']

    for sp in m.find('listofspecies').find_all('species'):
        cp = _reactome_id(sp, 'compartment')
        si = _reactome_id(sp, 'id')
        nm = sp.attrs['name']
        ids = []

        for i in sp.find('bqbiol:haspart').find_all('rdf:li'):
            ids.append(_reactome_res(i))

        ids = sorted(common.uniq_list(ids))
        species[si] = {'name': nm, 'comp': cp, 'ids': ids}

    for rea in m.find('listofreactions').find_all('reaction'):
        ri = _reactome_id(rea, 'id')
        refs = []

        for r in rea.find('bqbiol:isdescribedby').find_all('rdf:li'):
            refs.append(_reactome_res(r))

        refs = sorted(common.uniq_list(refs))
        reas = []

        for r in rea.find('listofreactants').find_all('speciesreference'):
            reas.append(_reactome_id(r, 'species'))

        reas = sorted(common.uniq_list(reas))
        prds = []

        for p in rea.find('listofproducts').find_all('speciesreference'):
            prds.append(_reactome_id(p, 'species'))

        prds = sorted(common.uniq_list(prds))
        note = rea.find('notes').text
        reactions[ri] = {
            'refs': refs,
            'reas': reas,
            'prds': prds,
            'note': note
        }

    return compartments, species, reactions


def _reactome_reactions_et():
    sbmlPfx = '{http://www.sbml.org/sbml/level2/version4}'
    compStr = '%scompartment' % sbmlPfx
    reacStr = '%sreaction' % sbmlPfx
    specStr = '%sspecies' % sbmlPfx
    species = {}
    compartments = {}
    reactions = {}
    sbmlfile = reactome_sbml()
    ctx = etree.iterparse(sbmlfile, events = ('end', ))

    for ev, elem in ctx:
        if elem.tag == compStr:
            k, v = _reactome_compartment(elem)
            compartments[k] = v

        elif elem.tag == reacStr:
            k, v = _reactome_reaction(elem)
            reactions[k] = v

        elif elem.tag == specStr:
            k, v = _reactome_species(elem)
            species[k] = v

        elem.clear()

        while elem.getprevious() is not None:
            del elem.getparent()[0]

    return compartments, species, reactions


def _reactome_compartment(elem):
    ci = _reactome_extract_id(elem.get('id'))
    nm = elem.get('name')

    return ci, nm


def _reactome_species(elem):
    bqBiolPfx = '{http://biomodels.net/biology-qualifiers/}'
    rdfPfx = '{http://www.w3.org/1999/02/22-rdf-syntax-ns#}'
    hasPartStr = '%shasPart' % bqBiolPfx
    resStr = '%sresource' % rdfPfx
    si = _reactome_extract_id(elem.get('id'))
    cp = _reactome_extract_id(elem.get('compartment'))
    nm = elem.get('name')
    ids = sorted(
        common.uniq_list(_reactome_collect_resources(elem, hasPartStr)))

    return si, {'name': nm, 'comp': cp, 'ids': ids}


def _reactome_reaction(elem):
    bqBiolPfx = '{http://biomodels.net/biology-qualifiers/}'
    rdfPfx = '{http://www.w3.org/1999/02/22-rdf-syntax-ns#}'
    sbmlPfx = '{http://www.sbml.org/sbml/level2/version4}'
    specStr = 'species'
    spRefStr = '%sspeciesReference' % sbmlPfx
    isDescStr = '%sisDescribedBy' % bqBiolPfx
    resStr = '%sresource' % rdfPfx
    lofReaStr = '%slistOfReactants' % sbmlPfx
    lofPrdStr = '%slistOfProducts' % sbmlPfx
    ri = _reactome_extract_id(elem.get('id'))
    refs = _reactome_collect_resources(elem, isDescStr)
    reas = _reactome_collect_species(elem, lofReaStr)
    prds = _reactome_collect_species(elem, lofPrdStr)
    note = elem.find('note').text  # prefix?

    return ri, {'refs': refs, 'reas': reas, 'prds': prds, 'note': note}


def _reactome_collect_resources(elem, tag):
    rdfPfx = '{http://www.w3.org/1999/02/22-rdf-syntax-ns#}'
    resStr = '%sresource' % rdfPfx
    liStr = '%sli' % rdfPfx
    res = []

    for i in elem.find('.//%s' % tag).iterfind('.//%s' % liStr):
        res.append(_reactome_extract_res(i.get(resStr)))

    return res


def _reactome_collect_species(elem, tag):
    sbmlPfx = '{http://www.sbml.org/sbml/level2/version4}'
    spRefStr = '%sspeciesReference' % sbmlPfx
    specStr = 'species'
    res = []

    for sp in elem.find('.//%s' % tag).iterfind('.//%s' % spRefStr):
        res.apped(_reactome_extract_id(sp.get(specStr)))

    return res


def get_acsn_effects():
    """
    Processes ACSN data, returns list of effects.
    """

    negatives = set(['NEGATIVE_INFLUENCE', 'UNKNOWN_NEGATIVE_INFLUENCE'])
    positives = set(
        ['TRIGGER', 'POSITIVE_INFLUENCE', 'UNKNOWN_POSITIVE_INFLUENCE'])
    directed = set([
        'UNKNOWN_TRANSITION', 'INTERACTION_TYPE', 'KNOWN_TRANSITION_OMITTED',
        'INHIBITION', 'UNKNOWN_POSITIVE_INFLUENCE', 'PROTEIN_INTERACTION',
        'UNKNOWN_CATALYSIS', 'POSITIVE_INFLUENCE', 'STATE_TRANSITION',
        'TRANSLATION', 'UNKNOWN_NEGATIVE_INFLUENCE', 'NEGATIVE_INFLUENCE',
        'MODULATION', 'TRANSCRIPTION', 'COMPLEX_EXPANSION', 'TRIGGER',
        'CATALYSIS', 'PHYSICAL_STIMULATION', 'UNKNOWN_INHIBITION', 'TRANSPORT'
    ])

    data = acsn_interactions()

    effects = []

    for l in data:
        if len(l) == 4:
            eff = set(l[2].split(';'))

            if len(eff & negatives) > 0:
                effects.append([l[0], l[1], '-'])

            elif len(eff & positives) > 0:
                effects.append([l[0], l[1], '+'])

            elif len(eff & directed) > 0:
                effects.append([l[0], l[1], '*'])

    return effects


def get_reactions(types = None, sources = None):
    if type(types) is list:
        types = set(types)

    if type(sources) is list:
        sources = set(sources)
    cachefile = os.path.join(
        settings.get('cachedir'),
        'reaction_interactions_by_source.pickle'
    )

    if os.path.exists(cachefile):
        interactions = pickle.load(open(cachefile, 'rb'))

    else:
        import pypath.utils.pyreact as pyreact
        rea = pyreact.PyReact()
        rea.load_all()
        rea.expand_by_source()
        interactions = rea.interactions_by_source

    for i in interactions:
        if (sources is None or i[4] in sources) and \
                (types is None or len(i[2] & types)):
            yield [
                i[0], i[1],
                ';'.join(list(i[2] if types is None else i[2] & types)),
                str(int(i[3])), i[4], ';'.join(list(i[5]))
            ]
