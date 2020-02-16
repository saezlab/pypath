#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
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

_logger = session_mod.Logger(name = 'dataio')
_log = _logger._log
_console = _logger._console

import urllib

try:
    import urllib2
except ImportError:
    # this works seemless in Py3:
    import urllib.request
    urllib2 = urllib.request

import httplib2
try:
    import urlparse
except:
    # this works seemless in Py3:
    import urllib.parse
    urlparse = urllib.parse

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
import xml.etree.cElementTree as ET
from lxml import etree
import time
import copy
import struct
import json
import csv
import pycurl
import webbrowser
import requests
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

if 'long' not in __builtins__:
    long = int

if 'unicode' not in __builtins__:
    unicode = str

CURSOR_UP_ONE = '\x1b[1A'
ERASE_LINE = '\x1b[2K'

# UniProt ID with isoform e.g. O14754-1
reupi = re.compile(r'([\w]{6,10})(?:-([0-9]{1,2}))?')

#
# thanks for http://stackoverflow.com/a/3239248/854988
#


CellPhoneDBAnnotation = collections.namedtuple(
    'CellPhoneDBAnnotation',
    (
        'receptor',
        'receptor_class',
        'peripheral',
        'secreted',
        'secreted_class',
        'transmembrane',
        'integrin',
    )
)


def read_table(cols,
               fileObject = None,
               data = None,
               sep = '\t',
               sep2 = None,
               rem = [],
               hdr = None,
               encoding = 'ascii'):
    """
    Generic function to read data tables.

    fileObject : file-like
        Any file like object: file opened for read, or StringIO buffer
    cols : dict
        Dictionary of columns to read. Keys identifying fields are returned
        in the result. Values are column numbers.
    sepLevel1 : str
        Field separator of the file.
    sepLevel2 : dict
        Subfield separators and prefixes.
        E.g. {2: ',', 3: '|'}
    hdr : int
        Number of header lines. If None, no headers assumed.
    rem : list
        Strings to remove. For each line these elements will be replaced with ''.
    """
    if data is None:
        if hasattr(fileObject, 'readline'):
            fileObject.seek(0)
        if hdr:
            for h in xrange(0, hdr):
                _ = next(fileObject)
        data = fileObject
    else:
        data = [l.strip() for l in data.split('\n') if len(l) > 0][hdr:]
    res = []
    for l in data:
        if type(l) is bytes:
            l = l.decode(encoding)
        for r in rem:
            l = l.replace(r, '')
        l = [f.strip() for f in l.split(sep)]
        if len(l) > max(cols.values()):
            dic = {}
            for name, col in iteritems(cols):
                field = l[col].strip()
                if sep2 is not None:
                    field = [
                        sf.strip() for sf in field.split(sep2) if len(sf) > 0
                    ]
                dic[name] = field
            res.append(dic)
    if fileObject is not None:
        fileObject.close()
    return res


def all_uniprots(organism = 9606, swissprot = None):

    return uniprot_input.all_uniprots(organism, swissprot)


def get_pdb():

    c = curl.Curl(urls.urls['uniprot_pdb']['url'], silent = False)
    data = c.result
    if data is None:
        return None, None
    data = data.split('\n')
    u_pdb = {}
    pdb_u = {}
    pdb = None
    pdb_re = re.compile(r'[0-9A-Z]{4}')
    for l in data:
        l = re.split('[ ]{2,}',
                     re.sub('[ ]+,[ ]+', ',', re.sub(r'[ ]*\(', '(', l)))
        if len(l[0]) == 4 and pdb_re.match(l[0]):
            pdb = l[0].lower()
            res = None if l[2] == '-' else float(l[2].replace(' A', ''))
            met = l[1]
        if pdb is not None and len(l) > 1:
            uniprots = l[1] if len(l) < 4 else l[3]
            uniprots = [
                u.split('(')[1].replace(')', '') for u in uniprots.split(',')
                if '(' in u
            ]
            pdb_u[pdb] = uniprots
            for u in uniprots:
                if u not in u_pdb:
                    u_pdb[u] = []
                u_pdb[u].append((pdb, met, res))

    return u_pdb, pdb_u


def pdb_complexes(organism = None):

    complexes = {}

    uniprot_pdb, pdb_uniprot = get_pdb_chains()
    del uniprot_pdb

    for pdb_id, chains in iteritems(pdb_uniprot):

        uniprots = tuple(chain['uniprot'] for chain in chains.values())

        if len(uniprots) == 1:

            continue

        # if the organism set and any of the UniProt IDs does not
        # belong to this organism we drop the complex
        if organism and reflists.is_not(uniprots, 'uniprot', organism):

            continue

        cplex = intera.Complex(
            components = uniprots,
            sources = 'PDB',
            ids = pdb_id,
        )

        if cplex.__str__() in complexes:

            complexes[cplex.__str__()] += cplex

        else:

            complexes[cplex.__str__()] = cplex

    return complexes


def corum_complexes(organism = 9606):

    annots = (
        'mithocondr',
        'nucleus',
        'endoplasmic reticulum',
        'cytoplasm',
        'transcriptional control',
        'vesicle docking',
        'extracellular matrix component',
        'cell-matrix adhesion',
        'cytokines',
        'cell death',
        'integrin receptor signalling pathway',
        'eukaryotic plasma membrane',
        'nuclear membrane',
        'cellular export and secretion',
        'cell-substrate adherens junction',
        'cytoskeleton',
        'receptor binding',
        'nucleolus',
        'transmembrane signal transduction',
        'transcription',
        'modification by phosphorylation',
        'cell-cell adhesion',
        'intercellular junction',
        'ion transport',
        'cell adhesion',
        'cell junction',
        'endocytosis',
    )

    organism = taxonomy.ensure_ncbi_tax_id(organism)

    complexes = {}

    c = curl.Curl(
        urls.urls['corum']['url'],
        silent = False,
        large = True,
        files_needed = ['allComplexes.txt'],
    )

    tab = csv.DictReader(c.result['allComplexes.txt'], delimiter = '\t')

    for rec in tab:

        cplex_organism = rec['Organism']

        if taxonomy.ensure_ncbi_tax_id(cplex_organism) != organism:

            continue

        uniprots = rec['subunits(UniProt IDs)'].split(';')

        pubmeds  = rec['PubMed ID'].split(';')
        name     = rec['ComplexName']

        cplex = intera.Complex(
            name = name,
            components = uniprots,
            sources = 'CORUM',
            references = pubmeds,
            ids = rec['ComplexID'],
            attrs = {
                'funcat': set(rec['FunCat description'].split(';')),
                'go': set(rec['GO description'].split(';')),
            },
        )

        if cplex.__str__() in complexes:

            complexes[cplex.__str__()].references.update(set(pubmeds))

        else:

            complexes[cplex.__str__()] = cplex

    return complexes


def complexportal_complexes(organism = 9606, return_details = False):
    """
    Complex dataset from IntAct.
    See more:
    http://www.ebi.ac.uk/intact/complex/
    http://nar.oxfordjournals.org/content/early/2014/10/13/nar.gku975.full.pdf
    """

    spec = {9606: 'Homo_sapiens'}

    zipurl = '%s/%s.zip' % (
        urls.urls['complex_portal']['url'],
        spec[organism],
    )
    c = curl.Curl(zipurl, large = True, silent = False)
    files = c.result

    errors = []
    complexes = {}
    details = []
    name_key = 'complex recommended name'

    for xmlname, xml in iteritems(c.result):

        soup = bs4.BeautifulSoup(xml, 'html.parser')
        interactors_xml = soup.find_all('interactor')
        interactors = {}
        interactions = {}

        for i in interactors_xml:

            if i.find('primaryref').attrs['db'] == 'uniprotkb':

                interactors[i.attrs['id']] = i.find('primaryref').attrs['id']

        interactions_xml = soup.find_all('interaction')

        for i in interactions_xml:

            description = ''
            pubmeds = []
            fullname = ''
            names = {}
            pdbs = []
            uniprots = []
            ids = collections.defaultdict(set)

            for a in i.find_all('attribute'):

                if a.attrs['name'] == 'curated-complex':
                    description = a.text

            for sr in i.find_all('secondaryref'):

                if sr.attrs['db'] == 'pubmed':
                    pubmeds.append(sr.attrs['id'])

                if sr.attrs['db'] == 'wwpdb':
                    pdbs.append(sr.attrs['id'])

            for pr in i.find_all('primaryref'):

                if pr.attrs['db'] in {'wwpdb', 'rcsb pdb', 'pdbe'}:
                    pdbs.append(pr.attrs['id'])

            for sr in i.find('xref').find_all('secondaryref'):

                if (
                    'reftype' in sr.attrs and
                    sr.attrs['db'] in {'intact', 'reactome'} and
                    sr.attrs['reftype'] == 'identity'
                ):

                    ids[sr.attrs['db']].add(sr.attrs['id'])

            pubmeds = list(set(pubmeds))
            pdbs = list(set(pdbs))
            fullname = (
                None
                    if i.find('fullname') is None else
                i.find('fullname').text
            )

            for a in i.find_all('alias'):

                names[a.attrs['type']] = a.text

            for intref in i.find_all('interactorref'):

                int_id = intref.text

                if int_id in interactors:

                    uniprot = interactors[int_id]

                    if uniprot.startswith('PRO'):

                        continue

                    uniprot = uniprot.split('-')[0]

                    uniprots.append(uniprot)

            if uniprots:

                if pdbs:

                    ids['PDB'].update(set(pdbs))

                cplex = intera.Complex(
                    components = uniprots,
                    name = names[name_key] if name_key in names else None,
                    references = set(pubmeds),
                    sources = 'ComplexPortal',
                    ids = ids,
                )

                if cplex.__str__() in complexes:

                    complexes[cplex.__str__()] += cplex

                else:

                    complexes[cplex.__str__()] = cplex

            details.append({
                'uniprots': uniprots,
                'pdbs': pdbs,
                'pubmeds': pubmeds,
                'fullname': fullname,
                'names': names,
                'description': description
            })

    if return_details:

        return complexes, details

    else:

        return complexes


def get_havugimana():
    """
    Downloads data from
    Supplement Table S3/1 from Havugimana 2012
    Cell. 150(5): 1068–1081.
    """

    url = urls.urls['havugimana']['url']
    c = curl.Curl(url, silent = False, large = True)
    fname = c.fileobj.name
    del c
    table = inputs_common.read_xls(fname)

    return table[3:]


def havugimana_complexes():
    """
    Retrieves complexes from
    Supplement Table S3/1 from Havugimana 2012
    Cell. 150(5): 1068–1081.
    """

    complexes = {}

    for rec in get_havugimana():

        cplex = intera.Complex(
            components = rec[2].split(','),
            sources = 'Havugimana2012',
            ids = rec[0],
        )

        complexes[cplex.__str__()] = cplex

    return complexes


def compleat_complexes(predicted = True):
    """
    Retrieves complexes from the Compleat database.
    """

    url = urls.urls['compleat']['url']
    c = curl.Curl(url, large = True, silent = False)
    tab = list(csv.DictReader(
        c.result,
        delimiter = '\t',
        fieldnames = (
            'compleat_id',
            'member_count',
            'predicted',
            'functions',
            'functions2',
            'nothing',
            'sources',
            'name',
            'method',
            'organisms',
            'pubmeds',
            'members',
        )
    ))

    complexes = {}

    for rec in tab:


        is_predicted = (
            rec['predicted'] and
            rec['predicted'].strip() == 'Predicted'
        )

        if is_predicted and not predicted:

            continue

        if not rec['members']:

            continue

        uniprots = []

        for entrez in rec['members'].split():

            uniprot = mapping.map_name0(entrez.strip(), 'entrez', 'uniprot')

            if uniprot:

                uniprots.append(uniprot)

        if not uniprots:

            continue

        name = rec['name']
        references = rec['pubmeds'].split(',') if rec['pubmeds'] else None
        sources = set(rec['sources'].split(',')) if is_predicted else set()
        sources.add('Compleat')

        cplex = intera.Complex(
            components = uniprots,
            sources = sources,
            references = references,
            name = name,
            ids = {'Compleat': rec['compleat_id']},
        )

        if cplex.__str__() in complexes:

            complexes[cplex.__str__()] += cplex

        else:

            complexes[cplex.__str__()] = cplex

    return complexes


def humap_complexes():

    url = urls.urls['proteincomplexes']['url']
    c = curl.Curl(url, large = True)

    complexes = {}

    for l in c.result:

        l = l.strip().split()

        for uniprots in itertools.product(*(
            mapping.map_name(entrez, 'entrez', 'uniprot') for entrez in l
        )):

            cplex = intera.Complex(
                components = uniprots,
                sources = 'hu.MAP',
            )

            complexes[cplex.__str__()] = cplex

    return complexes


def get_pdb_chains():

    def to_int(i):

        if i == 'None':

            return None

        return int(non_digit.sub('', i))

    c = curl.Curl(urls.urls['pdb_chains']['url'], silent = False)
    chains = c.result
    if chains is None:
        return None, None
    chains = chains.replace('\r', '').split('\n')
    del chains[0]
    del chains[0]
    pdb_u = {}
    u_pdb = {}
    non_digit = re.compile(r'[^\d.-]+')

    for l in chains:
        l = l.split('\t')
        if len(l) > 8:
            if l[0] not in pdb_u:
                pdb_u[l[0]] = {}
            pdb_u[l[0]][l[1]] = {
                'uniprot': l[2],
                'chain_beg': to_int(l[3]),
                'chain_end': to_int(l[4]),
                'pdb_beg': to_int(l[5]),
                'pdb_end': to_int(l[6]),
                'uniprot_beg': to_int(l[7]),
                'uniprot_end': to_int(l[8])
            }
            if (
                pdb_u[l[0]][l[1]]['pdb_end'] is not None and
                pdb_u[l[0]][l[1]]['pdb_beg'] is not None and
                pdb_u[l[0]][l[1]]['uniprot_beg'] is not None and
                pdb_u[l[0]][l[1]]['uniprot_end'] is not None and
                pdb_u[l[0]][l[1]]['pdb_end'] - pdb_u[l[0]][l[1]]['pdb_beg'] == \
                    pdb_u[l[0]][l[1]]['uniprot_end'] - pdb_u[l[0]][l[1]]['uniprot_beg']
            ):

                pdb_u[l[0]][l[1]]['offset'] = (pdb_u[l[0]][l[1]]['uniprot_beg']
                                               - pdb_u[l[0]][l[1]]['pdb_beg'])

            else:
                pdb_u[l[0]][l[1]]['offset'] = None

            if l[2] not in u_pdb:

                u_pdb[l[2]] = []

            u_pdb[l[2]].append({
                'pdb': l[0],
                'chain': l[1],
                'chain_beg': to_int(l[3]),
                'chain_end': to_int(l[4]),
                'pdb_beg': to_int(l[5]),
                'pdb_end': to_int(l[6]),
                'uniprot_beg': to_int(l[7]),
                'uniprot_end': to_int(l[8]),
                'offset': pdb_u[l[0]][l[1]]['offset']
            })

    return u_pdb, pdb_u


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


def domino_interactions():

    domino = get_domino()
    inter = []

    for l in domino:

        if (
            l[0] and
            l[1] and
            ''.join(l[5]) and
            ''.join([
                l[i]
                for i in range(10, 12) + range(14, 22) + range(24, 26)
            ]) and
            l[28] != '1'
        ):
            inter.append(l)

    return inter


def get_domino_ddi():

    domi = get_domino_ptms()
    return domi['ddi']


def get_domino_ptms():
    """
    The table comes from dataio.get_domino(), having the following fields:
    header = ['uniprot-A', 'uniprot-B', 'isoform-A', 'isoform-B', #3
    'exp. method', 'references', 'taxon-A', 'taxon-B', #7
    'role-A', 'role-B', 'binding-site-range-A', 'binding-site-range-B', #11
    'domains-A', 'domains-B', 'ptm-residue-A', 'ptm-residue-B', #15
    'ptm-type-mi-A', 'ptm-type-mi-B', 'ptm-type-A', 'ptm-type-B', #19
    'ptm-res-name-A', 'ptm-res-name-B', 'mutations-A', 'mutations-B', #23
    'mutation-effects-A', 'mutation-effects-B', 'domains-interpro-A', #26
    'domains-interpro-B', 'negative'] #28
    """
    domino = get_domino()
    try:
        miont = get_ontology('MI')
    except:
        miont = {}
    dmi = []
    ddi = []
    prg = progress.Progress(len(domino), 'Processing DOMINO', 11)
    for l in domino:
        prg.step()
        if (l[14].strip() != '' or l[15].strip() != '' or
            (l[10] != '' and l[11] != '')) and len(l[0]) > 0 and len(l[1]) > 0:
            uniprot1 = l[0]
            uniprot2 = l[1]
            # ptms
            if '-' not in l[14] and '-' not in l[15]:
                ptmre12 = [] if len(l[14]) == 0 else \
                    [int(x) for x in l[14].split(';')]
                ptmre21 = [] if len(l[15]) == 0 else \
                    [int(x) for x in l[15].split(';')]
                ptmty12 = [None for _ in ptmre12] if len(l[16]) == 0 else \
                    l[16].split(';')
                ptmty12 = [
                    None if x not in miont else miont[x] for x in ptmty12
                ]
                ptmrn12 = [None for _ in ptmre12] if len(l[20]) == 0 else \
                    l[20].split(';')
                ptmrn12 = [
                    None if x is None or x == '' or
                    len(x) < min(ptmre12[i] - 1, 11) else x[10]
                    if ptmre12[i] > 10 else x[ptmre12[i] - 1]
                    for i, x in enumerate(ptmrn12)
                ]
                ptmty21 = [None for _ in ptmre21] if len(l[17]) == 0 else \
                    l[17].split(';')
                ptmty21 = [
                    None if x not in miont else miont[x] for x in ptmty21
                ]
                ptmrn21 = [None for _ in ptmre21] if len(l[21]) == 0 else \
                    l[21].split(';')
                ptmrn21 = [
                    None if x is None or x == '' or
                    len(x) < min(ptmre21[i] - 1, 11) else x[10]
                    if ptmre21[i] > 10 else x[ptmre21[i] - 1]
                    for i, x in enumerate(ptmrn21)
                ]
                for i, resnum in enumerate(ptmre12):
                    res = intera.Residue(resnum, ptmrn12[i], uniprot2)
                    ptm = intera.Ptm(uniprot2,
                                     typ = ptmty12[i],
                                     residue = res,
                                     source = 'DOMINO')
                    dom = intera.Domain(uniprot1)
                    dm = intera.DomainMotif(
                        domain = dom,
                        ptm = ptm,
                        sources = 'DOMINO',
                        refs = l[5].split(';'))
            # binding sites
            if l[10] != '' and l[11] != '':
                try:
                    bssrt1 = [
                        int(x.split('-')[0]) for x in l[10].split(';')
                        if x != '' and x != '0'
                    ]
                    bsend1 = [
                        int(x.split('-')[1]) for x in l[10].split(';')
                        if x != '' and x != '0'
                    ]
                    bssrt2 = [
                        int(x.split('-')[0]) for x in l[11].split(';')
                        if x != '' and x != '0'
                    ]
                    bsend2 = [
                        int(x.split('-')[1]) for x in l[11].split(';')
                        if x != '' and x != '0'
                    ]
                except:
                    sys.stdout.write('Error processing line:\n')
                    sys.stdout.write(l)
                    sys.stdout.write('\n')
                    sys.stdout.flush()
                    return None
                bs1 = []
                bs2 = []
                if l[26] != '':
                    for i, n in enumerate(bssrt1):
                        bs1.append(
                            intera.Domain(
                                protein = uniprot1,
                                domain = l[26],
                                start = bssrt1[i],
                                end = bsend1[i],
                                domain_id_type = 'interpro',
                                isoform = l[2]))
                else:
                    for i, n in enumerate(bssrt1):
                        mot = intera.Motif(
                            protein = uniprot1,
                            start = bssrt1[i],
                            end = bsend1[i],
                            isoform = l[2])
                        bs1.append(
                            intera.Ptm(protein = uniprot1,
                                       motif = mot,
                                       source = 'DOMINO',
                                       isoform = l[2]))
                if l[27] != '':
                    for i, n in enumerate(bssrt2):
                        bs2.append(
                            intera.Domain(
                                protein = uniprot2,
                                domain = l[27],
                                start = bssrt2[i],
                                end = bsend2[i],
                                domain_id_type = 'interpro',
                                isoform = l[3]))
                else:
                    for i, n in enumerate(bssrt2):
                        mot = intera.Motif(
                            protein = uniprot2,
                            start = bssrt2[i],
                            end = bsend2[i],
                            isoform = l[3])
                        bs2.append(
                            intera.Ptm(
                                protein = uniprot2, motif = mot, source = 'DOMINO'))
                for one in bs1:
                    for two in bs2:
                        if one.__class__.__name__ == 'Domain' and \
                                two.__class__.__name__ == 'Domain':
                            dd = intera.DomainDomain(
                                one, two, sources = 'DOMINO')
                            ddi.append(dd)
                        if one.__class__.__name__ == 'Domain' and \
                                two.__class__.__name__ == 'Ptm':
                            dm = intera.DomainMotif(
                                domain = one,
                                ptm = two,
                                sources = 'DOMINO',
                                refs = l[6].split(';'))
                            dmi.append(dm)
                        if two.__class__.__name__ == 'Domain' and \
                                one.__class__.__name__ == 'Ptm':
                            dm = intera.DomainMotif(
                                domain = two,
                                ptm = one,
                                sources = 'DOMINO',
                                refs = l[6].split(';'))
                            dmi.append(dm)
    prg.terminate()
    return {'ddi': ddi, 'dmi': dmi}


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
    table = read_table(cols = cols, fileObject = buff, sep2 = subf, hdr = 1)
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
    table = read_table(cols = cols, fileObject = buff, sep = ',', hdr = 1)
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


def get_ontology(ontology):
    """
    Downloads an ontology using the bioservices module.
    """
    ols = bioservices.WSDLService("OLS", urls.urls['ols']['url'])
    ont = dict((x.key, x.value)
               for x in ols.serv.getAllTermsFromOntology(ontology).item)
    return ont


def get_listof_ontologies():
    """
    Returns a list of available ontologies using the bioservices module.
    """
    ols = bioservices.WSDLService("OLS", urls.urls['ols']['url'])
    olist = dict((x.key, x.value) for x in ols.serv.getOntologyNames().item)
    return olist


def residue_pdb(pdb, chain, residue):
    url = urls.urls['pdbsws']['url']
    params = urllib.urlencode({
        'plain': 1,
        'qtype': 'pdb',
        'id': pdb,
        'chain': chain,
        'res': residue
    })
    data = urllib2.urlopen(url + "?%s" % params)
    result = {}
    for l in data:
        if not l.startswith('//'):
            l = [x.strip() for x in l.split(':')]
            result[l[0]] = l[1]
    return result


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
        mapper = {}
        soup = bs4.BeautifulSoup(data.read(), 'html.parser')
        for block in soup.find_all('block'):
            seg = block.find_all('segment')
            chain = seg[0]['intobjectid'].split('.')[1]
            uniprot = seg[1]['intobjectid']
            pdbstart = int(non_digit.sub('', seg[0]['start']))
            pdbend = int(non_digit.sub('', seg[0]['end']))
            uniprotstart = int(non_digit.sub('', seg[1]['start']))
            uniprotend = int(non_digit.sub('', seg[1]['end']))
            if chain not in mapper:
                mapper[chain] = {}
            mapper[chain][pdbend] = {
                'uniprot': uniprot,
                'pdbstart': pdbstart,
                'uniprotstart': uniprotstart,
                'uniprotend': uniprotend
            }
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
                    if seg['pdbstart'] <= resnum:
                        offset = seg['uniprotstart'] - seg['pdbstart']
                        residue = {
                            'resnum': resnum + offset,
                            'offset': offset,
                            'uniprot': seg['uniprot'],
                            'chain': chain
                        }
                        return residue
        return None

    def clean(self):
        """
        Removes cached mappings, freeing up memory.
        """
        self.mappers = {}


def comppi_interaction_locations(organism = 9606):
    """
    Downloads and preprocesses protein interaction and cellular compartment
    association data from the ComPPI database.
    This data provides scores for occurrence of protein-protein interactions
    in various compartments.
    """

    ComppiLocation = collections.namedtuple(
        'ComppiLocation',
        [
            'location',
            'score',
        ],
    )

    ComppiInteraction = collections.namedtuple(
        'ComppiInteraction',
        [
            'id_a',
            'id_b',
            'loc_a',
            'loc_b',
        ],
    )

    def process_locations(loc):

        return tuple(
            ComppiLocation(location = llloc[0], score = float(llloc[1]))
            for llloc in
            (lloc.split(':') for lloc in loc.split('|'))
        )

    url = urls.urls['comppi']['url']
    post = {
        'fDlSet': 'comp',
        'fDlSpec': '0',
        'fDlMLoc': 'all',
        'fDlSubmit': 'Download'
    }
    c = curl.Curl(
        url,
        post = post,
        large = True,
        silent = False,
        compr = 'gz',
    )

    _ = next(c.result)

    for l in c.result:

        l = l.strip('\r\n').split('\t')

        organism_a = int(l[7])
        organism_b = int(l[15])

        if organism and (organism_a != organism or organism_b != organism):

            continue
        
        for uniprot1, uniprot2 in itertools.product(
            mapping.map_name(l[0], 'uniprot', 'uniprot'),
            mapping.map_name(l[8], 'uniprot', 'uniprot'),
        ):
            
            yield ComppiInteraction(
                id_a = uniprot1,
                id_b = uniprot2,
                loc_a = process_locations(l[2]),
                loc_b = process_locations(l[10]),
            )


def comppi_locations(organism = 9606, score_threshold = .7):

    result = collections.defaultdict(set)

    for iloc in comppi_interaction_locations(organism = organism):

        for label in ('a', 'b'):

            for loc in getattr(iloc, 'loc_%s' % label):

                if loc.location == 'N/A' or loc.score < score_threshold:

                    continue

                result[getattr(iloc, 'id_%s' % label)].add(loc)

    return dict(result)


def ramilowski_locations(long_notes = False):

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

    return result


def get_psite_phos(raw = True, organism = 'human', strict = True):
    """
    Downloads and preprocesses phosphorylation site data from PhosphoSitePlus.
    """

    url = urls.urls['psite_kin']['url']
    c = curl.Curl(
        url,
        silent = False,
        compr = 'gz',
        encoding = 'iso-8859-1',
        large = True,
    )
    orto = {}
    data = c.result
    cols = {
        'kinase': 2,
        'kinase_org': 3,
        'substrate': 6,
        'substrate_org': 8,
        'residue': 9,
        'motif': 11
    }
    data = read_table(cols = cols, fileObject = data, sep = '\t', hdr = 4)
    result = []
    non_digit = re.compile(r'[^\d.-]+')
    motre = re.compile(r'(_*)([A-Za-z]+)(_*)')
    
    for r in data:

        if organism is None or \
            ((r['kinase_org'] == organism or not strict) and \
            r['substrate_org'] == organism):

            if r['kinase_org'] != organism:
                korg = r['kinase_org']
                # attempting to map by orthology:
                if korg in taxonomy.taxa and organism in taxonomy.taxa:

                    ktaxid = taxonomy.taxa[korg]
                    taxid = taxonomy.taxa[organism]

                    if korg not in orto:
                        
                        orto[korg] = homology_mod.homologene_dict(
                            ktaxid,
                            taxid,
                            'refseqp',
                        )

                    korg_refseq = mapping.map_name(r['kinase'],
                                                    'uniprot',
                                                    'refseqp',
                                                    ktaxid)

                    kin_uniprot = \
                        list(
                            itertools.chain(
                                *map(
                                    lambda ors:
                                        mapping.map_name(ors,
                                                        'refseqp',
                                                        'uniprot',
                                                        taxid),
                                    itertools.chain(
                                        *map(
                                            lambda rs:
                                                orto[korg][rs],
                                            filter(
                                                lambda rs:
                                                    rs in orto[korg],
                                                korg_refseq
                                            )
                                        )
                                    )
                                )
                            )
                        )
            else:
                kin_uniprot = [r['kinase']]

            for kinase in kin_uniprot:

                r['resaa'] = r['residue'][0]
                r['resnum'] = int(non_digit.sub('', r['residue'][1:]))
                mot = motre.match(r['motif'])

                r['substrate'] = r['substrate'].split('_')[0] # excluding e.g. Q12809_VAR_014388
                sisoform = 1 if '-' not in r['substrate'] else \
                    int(r['substrate'].split('-')[1])
                r['substrate'] = r['substrate'].split('-')[0]

                kisoform = 1 if '-' not in kinase else int(kinase.split('-')[1])
                kinase = kinase.split('-')[0]

                r['substrate'] = r['substrate'].split('-')[0]

                if mot:
                    r['start'] = r['resnum'] - 7 + len(mot.groups()[0])
                    r['end'] = r['resnum'] + 7 - len(mot.groups()[2])
                    r['instance'] = r['motif'].replace('_', '').upper()
                else:
                    r['start'] = None
                    r['end'] = None
                    r['instance'] = None

                if raw:
                    r['kinase'] = kinase
                    result.append(r)
                else:
                    res = intera.Residue(r['resnum'], r['resaa'],
                                         r['substrate'],
                                         isoform = sisoform)

                    mot = intera.Motif(
                        r['substrate'],
                        r['start'],
                        r['end'],
                        instance = r['instance'],
                        isoform = sisoform)

                    ptm = intera.Ptm(protein = r['substrate'],
                                    residue = res,
                                    motif = mot,
                                    typ = 'phosphorylation',
                                    source = 'PhosphoSite',
                                    isoform = sisoform)

                    dom = intera.Domain(protein = kinase, isoform = kisoform)

                    dommot = intera.DomainMotif(
                        domain = dom, ptm = ptm, sources = ['PhosphoSite'])

                    result.append(dommot)

    return result

def ptm_orthology():
    """
    Returns an orthology translation dict of phosphosites
    based on phosphorylation sites table from PhosphoSitePlus.
    In the result all PTMs represented by a tuple of the following
    6 elements: UniProt ID, isoform (int), residue one letter code,
    residue number (int), NCBI Taxonomy ID (int), modification type.

    :param int source: Source taxon (NCBI Taxonomy).
    :param int target: Target taxon (NCBI Taxonomy).
    """
    result = {}

    nondigit = re.compile(r'[^\d]+')

    unknown_taxa = set([])

    for typ in common.psite_mod_types:

        groups = {}

        url = urls.urls['psite_%s' % typ[0]]['url']
        c = curl.Curl(url, silent = False, large = True)

        data = c.result

        for _ in xrange(4):

            __ = next(data)

        for r in data:

            r = r.decode('utf-8').split('\t')

            if len(r) < 10:

                continue

            uniprot = r[2]
            isoform = 1 if '-' not in uniprot else int(uniprot.split('-')[1])
            uniprot = uniprot.split('-')[0]
            aa = r[4][0]
            num = int(nondigit.sub('', r[4]))
            if r[6] not in taxonomy.taxa:
                unknown_taxa.add(r[6])
                continue

            tax = taxonomy.taxa[r[6]]
            group = int(r[5])

            this_site = (uniprot, isoform, aa, num, tax, typ[1])

            if group not in groups:
                groups[group] = set([])

            groups[group].add(this_site)

        for group, sites in iteritems(groups):

            for site1 in sites:

                for site2 in sites:

                    if site1[4] == site2[4]:

                        continue

                    if site1 not in result:

                        result[site1] = {}

                    if site2[4] not in result[site1]:

                        result[site1][site2[4]] = set([])

                    result[site1][site2[4]].add(site2)

    if len(unknown_taxa):
        sys.stdout.write('\t:: Unknown taxa encountered:\n\t   %s\n' %
                         ', '.join(sorted(unknown_taxa)))

    return result

def get_psite_p(organism = 'human'):
    """
    Downloads the phosphorylation site dataset from PhosphoSitePlus.
    """
    result = []
    url = urls.urls['psite_p']['url']
    nondigit = re.compile(r'[^\d]+')
    remot = re.compile(r'(_*)([A-Za-z]+)(_*)')

    c = curl.Curl(url, silent = False, large = True)
    data = c.result

    for _ in xrange(4):
        null = c.result.readline()

    for r in data:

        r = r.split('\t')

        if len(r) > 9 and (organism is None or r[6] == organism):

            uniprot = r[2]
            isoform = 1 if '-' not in uniprot else int(uniprot.split('-')[1])
            uniprot = uniprot.split('-')[0]
            typ = r[3].lower()
            if len(typ) == 0:
                typ = r[4].split('-')[1] if '-' in r[4] else None
            aa = r[4][0]
            num = int(nondigit.sub('', r[4]))
            motif = remot.match(r[9])
            if motif:
                start = num - 7 + len(motif.groups()[0])
                end = num + 7 - len(motif.groups()[2])
                instance = r[9].replace('_', '').upper()
            else:
                start = None
                end = None
                instance = None

            res = intera.Residue(num, aa, uniprot, isoform = isoform)
            mot = intera.Motif(
                uniprot, start, end, instance = instance, isoform = isoform)
            ptm = intera.Ptm(uniprot,
                             typ = typ,
                             motif = mot,
                             residue = res,
                             source = 'PhosphoSite',
                             isoform = isoform)
            result.append(ptm)

    return result

def get_psite_reg():
    """
    Downloads and preprocesses the regulatory sites dataset from
    PhosphoSitePlus. This data provides information about which
    proteins a PTM disrupts or induces the interaction with.
    """
    kwds_pos = {
        'enzymatic activity, induced',
        'activity, induced',
        'protein stabilization',
        'receptor inactivation, inhibited',
        'receptor desensitization, inhibited',
        'receptor internalization, inhibited',
        'receptor recycling, induced'
    }

    kwds_neg = {
        'enzymatic activity, inhibited',
        'activity, inhibited',
        'protein degradation',
        'receptor inactivation, induced',
        'receptor desensitization, induced',
        'receptor internalization, induced',
        'receptor recycling, inhibited'
    }

    url = urls.urls['psite_reg']['url']
    c = curl.Curl(url, silent = False, compr = 'gz',
                  encoding = 'iso-8859-1', large = True)
    data = c.result
    cols = {
        'uniprot': 3,
        'organism': 6,
        'mod': 7,
        'on_function': 11,
        'on_process': 12,
        'on_interact': 13,
        'pmids': 15,
        'comments': 19
    }

    data = read_table(cols = cols, fileObject = data, sep = '\t', hdr = 4)
    regsites = {}

    for r in data:

        interact = [[y.replace(')', '').strip() for y in x.split('(')]
                    for x in r['on_interact'].strip().split(';') if len(x) > 0]
        induces = [x[0] for x in interact if x[1] == 'INDUCES']
        disrupts = [x[0] for x in interact if x[1] == 'DISRUPTS']
        mod = r['mod']
        modt = r['mod'].split('-')
        mod = list(modt[0])
        aa = mod.pop(0)
        modt = modt[1]
        res = ''.join(mod)
        isoform = (
            int(r['uniprot'].split('-')[1])
            if '-' in r['uniprot']
            else 1
        )
        uniprot = r['uniprot'].split('-')[0]

        if uniprot not in regsites:
            regsites[uniprot] = []

        function = set(map(lambda f: f.strip(),
                           r['on_function'].split(';')))

        regsites[uniprot].append({
            'aa':       aa,
            'res':      res,
            'modt':     modt,
            'organism': r['organism'],
            'pmids':    set(map(lambda f: f.strip(),
                                r['pmids'].split(';'))),
            'induces':  induces,
            'disrupts': disrupts,
            'isoform':  isoform,
            'function': function,
            'process':  set(map(lambda f: f.strip(),
                                r['on_process'].split(';'))),
            'comments': r['comments'],
            'positive': bool(kwds_pos & function),
            'negative': bool(kwds_neg & function)
        })

    return regsites

def regsites_one_organism(organism = 9606):
    """
    Returns PhosphoSitePlus regulatory sites translated to
    one organism by orthology. Residue numbers will be translated
    where necessary, while gene symbols will be translated to
    UniProt IDs of the given organism.
    This works with human, mouse or rat.

    :param int organism:
        NCBI Taxonomy ID of the target organism. In this
        method possible values are human, mouse or rat, as these species
        provide the vast majority of the data, and are close enough to each
        other that the sites can be safely translated between orthologous
        proteins by sequence alignement.
    """

    def genesymbols2uniprots(genesymbols, tax):
        return (
            set(
                itertools.chain(
                    *map(
                        lambda gs:
                            mapping.map_name(
                                gs,
                                'genesymbol',
                                'uniprot',
                                ncbi_tax_id = tax,
                            ),
                        genesymbols
                    )
                )
            )
        )

    def translate_uniprots(uniprots, homo):
        return (
            set(
                itertools.chain(
                    *map(
                        lambda usrc:
                            homo[usrc] if usrc in homo else [],
                        uniprots
                    )
                )
            )
        )

    result = {}

    organisms = set([9606, 10090, 10116])

    mod_types = dict(common.psite_mod_types2)

    regsites = get_psite_reg()

    other_organisms = organisms - set([organism])

    homology = (
        dict(
            map(
                lambda other:
                    (
                        other,
                        homology_mod.homologene_uniprot_dict(
                            source = other,
                            target = organism,
                        )
                    ),
                other_organisms
            )
        )
    )

    ptm_homology = ptm_orthology()

    proteome = uniprot_input.all_uniprots(organism = organism, swissprot = 'YES')

    for substrate, regs in iteritems(regsites):

        subs = []

        if substrate in proteome:
            subs = [substrate]
        else:
            for other, homo in iteritems(homology):
                if substrate in homo:
                    subs = homo[substrate]

        for sub in subs:

            if sub not in result:
                result[sub] = {}

            for reg in regs:

                reg_organism = taxonomy.taxa[reg['organism']]

                if reg_organism not in organisms:
                    continue

                mod_type = mod_types[reg['modt']]
                resnum = int(reg['res'])

                psite_key = (
                    substrate,
                    reg['isoform'],
                    reg['aa'],
                    resnum,
                    reg_organism,
                    mod_type,
                )

                if reg_organism != organism:

                    regs_target = []
                    disrupts    = []
                    induces     = []

                    if psite_key in ptm_homology:

                        if organism in ptm_homology[psite_key]:

                            regs_target = ptm_homology[psite_key][organism]

                    if len(regs_target):

                        disrupts = genesymbols2uniprots(
                            reg['disrupts'],
                            reg_organism,
                        )
                        disrupts = translate_uniprots(
                            disrupts,
                            homology[reg_organism],
                        )
                        induces  = genesymbols2uniprots(
                            reg['induces'],
                            reg_organism,
                        )
                        induces  = translate_uniprots(
                            induces,
                            homology[reg_organism],
                        )

                else:

                    regs_target = [psite_key]

                    disrupts = genesymbols2uniprots(reg['disrupts'], organism)
                    induces  = genesymbols2uniprots(reg['induces'], organism)

                for regt in regs_target:

                    modkey = (regt[2], regt[3], regt[5])

                    if modkey not in result[sub]:

                        result[sub][modkey] = {
                            'induces':  set([]),
                            'disrupts': set([]),
                            'pmids':    set([]),
                            'isoforms': set([]),
                            'process':  set([]),
                            'function': set([]),
                            'positive': False,
                            'negative': False,
                            'comments': []
                        }

                    result[sub][modkey]['induces'].update(induces)
                    result[sub][modkey]['disrupts'].update(disrupts)
                    result[sub][modkey]['process'].update(reg['process'])
                    result[sub][modkey]['function'].update(reg['function'])
                    result[sub][modkey]['isoforms'].update([regt[1]])
                    result[sub][modkey]['pmids'].update(reg['pmids'])
                    result[sub][modkey]['positive'] = \
                        result[sub][modkey]['positive'] or reg['positive']
                    result[sub][modkey]['negative'] = \
                        result[sub][modkey]['negative'] or reg['negative']
                    if len(reg['comments']):
                        result[sub][modkey]['comments'].append(reg['comments'])


    return result

def regsites_tab(regsites, outfile = None):
    """
    Exports PhosphoSite regulatory sites as a tabular file, all
    IDs translated to UniProt.
    """
    header = [
        'uniprot_a', 'isoform_a', 'a_res_aa', 'a_res_num', 'a_mod_type',
        'effect', 'uniprot_b', 'references'
    ]
    result = []
    for uniprot, regsite in iteritems(regsites):
        isoform = '1'
        uniprot = uniprot.split('-')
        if len(uniprot) > 1:
            isoform = uniprot[1]
        uniprot = uniprot[0]
        for r in regsite:
            if r['organism'] == 'human':
                for i in r['induces']:
                    other = mapping.map_name(i, 'genesymbol', 'uniprot')
                    for o in other:
                        if o != 'unmapped':
                            result.append([
                                uniprot, isoform, r['aa'], r['res'], r['modt'],
                                '+', o
                            ])
                for i in r['disrupts']:
                    other = mapping.map_name(i, 'genesymbol', 'uniprot')
                    for o in other:
                        if o != 'unmapped':
                            result.append([
                                uniprot, isoform, r['aa'], r['res'], r['modt'],
                                '-', o, ';'.join(r['pmids'])
                            ])
    if outfile is not None:
        out = '\t'.join(header) + '\n'
        for r in result:
            out += '\t'.join(r) + '\n'
        with open(outfile, 'w') as f:
            f.write(out)
    return result


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


def get_pepcyber(cache = None):
    """
    Downloads phosphoprotein binding protein interactions
    from the PEPCyber database.
    """

    def get_cells(row):
        cells = row.find_all('td')
        if len(cells) == 10:
            sp = cells[4].find('span')
            if sp is not None and 'class' in sp.attrs \
                    and 'sequence' in sp.attrs['class']:
                return cells

    url = urls.urls['pepcyber']['url']
    # this is huge, takes a few minutes!
    c = curl.Curl(url, silent = False, timeout = 600, encoding = 'iso-8859-1')
    data = c.result
    soup = bs4.BeautifulSoup(data, 'html.parser')
    rows = soup.find_all('tr')
    result = []
    uniprots = {}
    if cache is None:
        cache = os.path.join(settings.get('cachedir'), 'pepcyber-uniprots')
    if os.path.exists(cache):
        with open(cache, 'r') as f:
            for l in f:
                l = l.split('\t')
                if l[0] == u'\xce':
                    continue
                l.extend(['', ''])
                uniprots[l[0].strip()] = [l[1].strip(), l[2].strip()]
    prg = progress.Progress(len(rows), 'Retrieving and processing data', 7)
    notfound = []

    for row in rows:
        prg.step()
        cells = get_cells(row)
        if cells is None:
            continue

        thisRow = [c.text.strip() for c in cells]
        if len(thisRow) > 9 and thisRow[5].isdigit():
            inum = int(row.find('a')['name'])
            thisRow[9] = None if 'p' not in thisRow[4] else \
                thisRow[4][thisRow[4].index('p') + 1]
            if thisRow[2] not in uniprots or thisRow[3] not in uniprots:
                up = pepcyber_uniprot(inum)
                uniprots.update(up)
            if thisRow[2] in uniprots and thisRow[3] in uniprots:
                thisRow.extend(uniprots[thisRow[2]])
                thisRow.extend(uniprots[thisRow[3]])
                result.append(thisRow[1:])
            else:
                notfound.append([thisRow[2], thisRow[3], inum])
    prg.terminate()
    with open(cache, 'w') as f:
        for g, u in iteritems(uniprots):
            if g[0] == u'\xce':
                continue
            f.write('\t'.join([g] + u) + '\n')
    return result


def pepcyber_uniprot(num):
    result = {}
    url = urls.urls['pepcyber_details']['url'] % num
    c = curl.Curl(url, cache = False, encoding = 'iso-8859-1')
    data = c.result
    if data is None:
        return result
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
            if len(gname) > 0 and gname[0] != u'\xce':
                result[gname] = [swprot, refseq]
            gname = None
        prev = td.text.strip()
    return result


def pdzbase_interactions():
    """
    Downloads data from PDZbase. Parses data from the HTML tables.
    """
    
    PDZbaseInteraction = collections.namedtuple(
        'PDZbaseInteraction',
        [
            'uniprot_pdz',
            'isoform_pdz',
            'uniprot_ligand',
            'isoform_ligand',
            'genesymbol_pdz',
            'genesymbol_ligand',
            'pdz_domain',
            'organism',
            'pubmed',
        ],
    )
    
    
    url = urls.urls['pdzbase']['url_rescued']
    c = curl.Curl(url, silent = False)
    data = c.result
    soup = bs4.BeautifulSoup(data, 'html.parser')
    rows = (
        soup.find_all('table')[3].find('table').find('table').find_all('tr')
    )
    result = []
    
    del rows[0]
    
    for r in rows:
        
        r = [c.text.strip() for c in r.find_all('td')]
        
        uniprot_pdz, isoform_pdz = reupi.match(r[1]).groups()
        uniprot_ligand, isoform_ligand = reupi.match(r[4]).groups()
        
        result.append(
            PDZbaseInteraction(
                uniprot_pdz = uniprot_pdz,
                isoform_pdz = int(isoform_pdz) if isoform_pdz else 1,
                uniprot_ligand = uniprot_ligand,
                isoform_ligand = int(isoform_ligand) if isoform_ligand else 1,
                genesymbol_pdz = r[0],
                genesymbol_ligand = r[3],
                pdz_domain = int(r[2]),
                organism = taxonomy.ensure_ncbi_tax_id(r[5]),
                pubmed = int(r[6]),
            )
        )
    
    return result


def get_domino(none_values = False, outfile = None):
    
    result = []
    taxid = re.compile(r'taxid:(.*)\([a-zA-Z ]*\)')
    miont = re.compile(r'MI:[0-9]{4}\((.*)\)')
    binds = re.compile(r'([-0-9]*);.*')
    domai = re.compile(r'.*;.*;.*\((.*)\)')
    dipro = re.compile(r'.*;.*;.+:(IPR[0-9]*).*')
    ptmrs = re.compile(r'([-0-9]*);.*')
    ptmmi = re.compile(r'[0-9]*;(MI:[0-9]*)\(.*\);.*;.*')
    ptmrn = re.compile(
        r'.*sequence:[\s]*[0-9]+-[0-9]+[\s]*:[\s]*([A-Z]{10,}).*')
    ptmty = re.compile(r'[0-9]*;MI:[0-9]*\((.*)\);.*;.*')
    refrs = re.compile(r'(pubmed|doi):["]*([-0-9a-zA-Z\.\(\)/]*)["]*')
    url = urls.urls['domino']['rescued']
    c = curl.Curl(url, silent = False, large = True)
    data = c.result
    _ = next(data)
    
    header = [
        'uniprot-A', 'uniprot-B', 'isoform-A', 'isoform-B', 'exp. method',
        'references', 'taxon-A', 'taxon-B', 'role-A', 'role-B',
        'binding-site-range-A', 'binding-site-range-B', 'domains-A',
        'domains-B', 'ptm-residue-A', 'ptm-residue-B', 'ptm-type-mi-A',
        'ptm-type-mi-B', 'ptm-type-A', 'ptm-type-B', 'ptm-res-name-A',
        'ptm-res-name-B', 'mutations-A', 'mutations-B', 'mutation-effects-A',
        'mutation-effects-B', 'domains-interpro-A', 'domains-interpro-B',
        'negative'
    ]
    
    for r in data:
        
        r = r.strip().split('\t')
        if len(r) < 39:
            continue
        thisRow = [
            None if ':' not in r[0] else r[0].split(':')[1].split('-')[0], None
            if ':' not in r[1] else r[1].split(':')[1].split('-')[0], '1'
            if '-' not in r[0] else r[0].split('-')[1], '1'
            if '-' not in r[1] else r[1].split('-')[1], None if
            miont.match(r[6]) is None else miont.match(r[6]).groups(1)[0], None
            if refrs.match(r[8]) is None else refrs.match(r[8]).groups(1)[1],
            None if taxid.match(r[9]) is None else
            taxid.match(r[9]).groups(1)[0], None if taxid.match(r[10]) is None
            else taxid.match(r[10]).groups(1)[0], None
            if miont.match(r[11]) is None else miont.match(r[11]).groups(1)[0],
            None if miont.match(r[16]) is None else
            miont.match(r[17]).groups(1)[0], ';'.join([
                '' if binds.match(x) is None else binds.match(x).groups(1)[0]
                for x in r[32].split(',')
            ]), ';'.join([
                '' if binds.match(x) is None else binds.match(x).groups(1)[0]
                for x in r[33].split(',')
            ]), ';'.join([
                '' if domai.match(x) is None else domai.match(x).groups(1)[0]
                for x in r[32].split(',')
            ]), ';'.join([
                '' if domai.match(x) is None else domai.match(x).groups(1)[0]
                for x in r[33].split(',')
            ]), ';'.join([
                '' if ptmrs.match(x) is None else ptmrs.match(x).groups(1)[0]
                for x in r[34].split('|')
            ]), ';'.join([
                '' if ptmrs.match(x) is None else ptmrs.match(x).groups(1)[0]
                for x in r[35].split('|')
            ]), ';'.join([
                '' if ptmmi.match(x) is None else ptmmi.match(x).groups(1)[0]
                for x in r[34].split('|')
            ]), ';'.join([
                '' if ptmmi.match(x) is None else ptmmi.match(x).groups(1)[0]
                for x in r[35].split('|')
            ]), ';'.join([
                '' if ptmty.match(x) is None else ptmty.match(x).groups(1)[0]
                for x in r[34].split('|')
            ]), ';'.join([
                '' if ptmty.match(x) is None else ptmty.match(x).groups(1)[0]
                for x in r[35].split('|')
            ]), ';'.join([
                '' if ptmrn.match(x) is None else ptmrn.match(x).groups(1)[0]
                for x in r[34].split('|')
            ]), ';'.join([
                '' if ptmrn.match(x) is None else ptmrn.match(x).groups(1)[0]
                for x in r[35].split('|')
            ]), ';'.join([
                '' if ptmrs.match(x) is None else ptmrs.match(x).groups(1)[0]
                for x in r[36].split('|')
            ]), ';'.join([
                '' if ptmrs.match(x) is None else ptmrs.match(x).groups(1)[0]
                for x in r[37].split('|')
            ]), ';'.join([
                '' if ptmty.match(x) is None else ptmty.match(x).groups(1)[0]
                for x in r[36].split('|')
            ]), ';'.join([
                '' if ptmty.match(x) is None else ptmty.match(x).groups(1)[0]
                for x in r[37].split('|')
            ]), '' if dipro.match(r[32]) is None else
            dipro.match(r[32]).groups(1)[0], '' if dipro.match(r[33]) is None
            else dipro.match(r[33]).groups(1)[0], '0'
            if r[38].strip() == '-' else '1'
        ]
        if not none_values:
            thisRow = ['' if x is None else x for x in thisRow]
        result.append(thisRow)
    if outfile:
        with open(outfile, 'w') as outf:
            outf.write('\t'.join(header) + '\n')
            for r in result:
                outf.write('\t'.join(['' if x is None else x
                                      for x in r]) + '\n')
    return result


def get_elm_domains():
    result = {}
    url = urls.urls['ielm_domains']['url']
    c = curl.Curl(url, silent = False)
    data = c.result
    soup = bs4.BeautifulSoup(data, 'html.parser')
    tbl = soup.find('table').find_all('td')
    rows = [tbl[x:x + 4] for x in xrange(0, len(tbl), 4)]
    for r in rows:
        uniprot = r[1].text
        motif = r[0].text
        if uniprot not in result:
            result[uniprot] = {}
        if motif not in result[uniprot]:
            result[uniprot][motif] = []
        result[uniprot][motif].append((r[2].text, r[3].text))
    return result


def get_phosphoelm(organism = 9606, ltp_only = True):
    """
    Downloads kinase-substrate interactions from phosphoELM.
    Returns list of dicts.

    :param int organism: NCBI Taxonomy ID.
    :param bool ltp_only: Include only low-throughput interactions.
    """
    result = []
    non_digit = re.compile(r'[^\d.-]+')

    if organism is None:
        _organism = None
    elif organism in taxonomy.phosphoelm_taxids:
        _organism = taxonomy.phosphoelm_taxids[organism]
    else:
        sys.stdout.write('\t:: Unknown organism: `%u`.\n' % organism)
        return []

    url = urls.urls['p_elm']['url']
    c = curl.Curl(url, silent = False)
    data = c.result
    data = [
        n for d, n in iteritems(data)
        if d.startswith(urls.urls['p_elm']['psites'])
    ]
    data = data[0] if len(data) > 0 else ''
    data = [l.split('\t') for l in data.split('\n')]
    kinases = get_phelm_kinases()
    del data[0]

    for l in data:

        if len(l) == 9 and (l[7] == _organism or _organism is None) \
            and (not ltp_only or l[6] == 'LTP'):

            l[1] = 1 if '-' not in l[0] else int(l[0].split('-')[1])
            l[0] = l[0].split('-')[0]
            del l[-1]

            if len(l[5]) > 0 and l[5] in kinases:
                kinase = kinases[l[5]]

                result.append({
                    'instance': None,
                    'isoform': l[1],
                    'resaa': l[3],
                    'resnum': int(non_digit.sub('', l[2])),
                    'start': None,
                    'end': None,
                    'substrate': l[0],
                    'kinase': kinase,
                    'references': l[4].split(';'),
                    'experiment': l[6],
                    'organism': l[7]
                })

    return result


def phelm_interactions(organism = 'Homo sapiens'):
    result = []
    data = get_phosphoelm(ltp_only = True)
    for l in data:
        result.append([
            l['kinase'], l['substrate'], ';'.join(l['references']),
            l['organism']
        ])
    return result


def phelm_psites():
    result = []
    data = get_phosphoelm()
    kinases = get_phelm_kinases()
    for l in data:
        l.append('1' if '-' not in l[0] else l[0].split('-')[1])
        l[0] = l[0].split('-')[0]
        l.append('' if l[4] not in kinases else kinases[l[4]])
        result.append(l)
    return result


def get_phelm_kinases():
    result = {}
    url = urls.urls['p_elm_kin']['url']
    c = curl.Curl(url, silent = False)
    data = c.result
    soup = bs4.BeautifulSoup(data, 'html.parser')
    for row in soup.find_all('table')[1].find_all('tr'):
        thisRow = [x.text for x in row.find_all('td')]
        if len(thisRow) > 2 and len(thisRow[2].strip()) > 0:
            result[thisRow[0]] = thisRow[2].strip()
    return result


def get_elm_classes():
    url = urls.urls['elm_class']['url']
    c = curl.Curl(url, silent = False)
    data = c.result
    data = [
        x.split('\t')[:-2] for x in data.replace('"', '').split('\n')[6:]
        if len(x) > 0
    ]
    return dict(zip([x[1] for x in data], data))


def get_elm_instances():
    
    url = urls.urls['elm_inst']['url']
    c = curl.Curl(url, silent = False)
    data = c.result
    data = data.replace('"', '').split('\t')
    data = data[6:]


def elm_interactions():
    """
    Downlods manually curated interactions from ELM.
    This is the gold standard set of ELM.
    """
    
    
    def number_or_none(value, typ = int):
        
        return typ(value) if value != 'None' else None
    
    
    retax = re.compile(r'"([0-9]+)"\([-:/,\.\[\]\(\)\w\s]+\)')
    
    ELMInteraction = collections.namedtuple(
        'ELMInteraction',
        [
            'motif_elm',
            'domain_pfam',
            'uniprot_motif',
            'uniprot_domain',
            'isoform_motif',
            'isoform_domain',
            'start_motif',
            'end_motif',
            'start_domain',
            'end_domain',
            'affinity_min',
            'affinity_max',
            'pubmeds',
            'taxon_motif',
            'taxon_domain',
        ],
    )
    
    result = []
    url = urls.urls['elm_int']['url']
    c = curl.Curl(url, silent = False)
    data = c.result
    data = data.split('\n')
    del data[0]
    
    for l in data:
        
        if not l:
            
            continue
        
        l = tuple(x.strip() for x in l.split('\t'))
        
        uniprot_mofif, isoform_motif = reupi.match(l[2]).groups()
        uniprot_domain, isoform_domain = reupi.match(l[3]).groups()
        
        result.append(
            ELMInteraction(
                motif_elm = l[0],
                domain_pfam = l[1],
                uniprot_motif = uniprot_mofif,
                uniprot_domain = uniprot_domain,
                isoform_motif = int(isoform_motif) if isoform_motif else 1,
                isoform_domain = int(isoform_domain) if isoform_domain else 1,
                start_motif = int(l[4]),
                end_motif = int(l[5]),
                start_domain = number_or_none(l[6]),
                end_domain = number_or_none(l[7]),
                affinity_min = number_or_none(l[8], float),
                affinity_max = number_or_none(l[9], float),
                pubmeds = tuple(map(int, l[10].split(','))) if l[10] else (),
                taxon_motif = int(retax.match(l[11]).groups()[0]),
                taxon_domain = int(retax.match(l[12]).groups()[0]),
            )
        )
    
    return result


def pfam_uniprot(uniprots, infile = None):
    
    result = {}
    url = urls.urls['pfam_up']['url']
    infile = infile if infile is not None \
        else os.path.join(settings.get('cachedir'), 'pfam-regions.tab.gz')
    if not os.path.exists(infile):
        sys.stdout.write('\t:: Downloading data...\n')
        sys.stdout.flush()
        if hasattr(urllib, 'urlretrieve'):
            urllib.urlretrieve(url, infile)
        else:
            _urllib.request.urlretrieve(url, infile)
    sys.stdout.write('\t:: Processing domains...\n')
    sys.stdout.flush()
    gzfile = gzip.open(infile, mode = 'r')
    prg = progress.Progress(len(uniprots), 'Looking up domains', 1)
    for l in gzfile:
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


def get_dbptm(organism = 9606):

    if organism is None:
        _organism = None
    elif organism in taxonomy.dbptm_taxids:
        _organism = taxonomy.dbptm_taxids[organism]
    else:
        sys.stdout.write('\t:: Unknown organism: `%u`.\n' % organism)
        return []
    
    url = urls.urls['dbptm']['old_table']
    c = curl.Curl(url, silent = False, large = True)
    data = []

    hdr = next(c.result).strip().split('\t')

    for l in c.result:

        l = l.strip().split('\t')

        data.append(dict(
            (
                key,
                (
                    None
                        if val == '' else
                    val.split(';')
                        if key in {'references', 'kinase'} else
                    int(val)
                        if val.isdigit() else
                    val
                )
            )
            for key, val in zip(hdr, l)
        ))

    return data


def get_dbptm_old(organism = 9606):
    """
    Downloads enzyme-substrate interactions from dbPTM.
    Returns list of dicts.
    """
    if organism is None:
        _organism = None
    elif organism in taxonomy.dbptm_taxids:
        _organism = taxonomy.dbptm_taxids[organism]
    else:
        sys.stdout.write('\t:: Unknown organism: `%u`.\n' % organism)
        return []

    result = []
    byre = re.compile(r'.*by\s([A-Za-z0-9\s]+)\.*')
    andre = re.compile(r',|and')
    non_digit = re.compile(r'[^\d.-]+')

    for url in urls.urls['dbptm']['urls']:

        c = curl.Curl(url, silent = False)
        extra = c.result

        for k, data in iteritems(extra):

            data = [x.split('\t') for x in data.split('\n')]

            for l in data:

                if len(l) > 8:

                    if _organism:
                        mnemonic = l[0].split('_')[1].strip()
                        if mnemonic != _organism:
                            continue

                    resnum = int(non_digit.sub('', l[2]))

                    ptm = ({
                        'substrate': l[1],
                        'typ': l[7].lower(),
                        'resaa': l[8][6],
                        'resnum': resnum,
                        'instance': l[8].strip(),
                        'references': l[4].split(';'),
                        'source': l[5].split()[0],
                        'kinase': None if byre.match(l[3]) is None else [
                            i.strip()
                            for i in andre.split(
                                byre.match(l[3]).groups(1)[0])
                        ],
                        'start': resnum - 6,
                        'end': resnum + 6
                    })

                    if ptm['kinase'] is not None:

                        if 'autocatalysis' in ptm['kinase']:

                            ptm['kinase'].append(ptm['substrate'])
                            ptm['kinase'].remove('autocatalysis')

                        ptm['kinase'] = [
                            k.replace('host', '').strip()
                            for k in ptm['kinase']
                        ]

                        ptm['kinase'] = [
                            k for k in ptm['kinase'] if len(k) > 0
                        ]

                        if len(ptm['kinase']) == 0:
                            ptm['kinase'] = None

                    result.append(ptm)

    return result


def dbptm_interactions():
    result = []
    data = get_dbptm()
    for r in data:
        if r['kinase'] is not None:
            for src in r['kinase']:
                result.append([
                    src, r['substrate'],
                    ';'.join([i for i in r['references'] if i != '-'])
                ])
    return result


def get_phosphonetworks():
    result = []
    reres = re.compile(r'([A-Z])([0-9]+)')
    non_digit = re.compile(r'[^\d.-]+')
    motre = re.compile(r'(-*)([A-Za-z]+)(-*)')
    url = urls.urls['phosnw']['url']
    c = curl.Curl(url, silent = False)
    data = c.result
    if data is None:
        return None
    data = data.split('\n')
    for l in data:
        if l.startswith('>'):
            substrate = l[1:].strip()
        elif len(l.split('\t')) >= 4:
            l = [x.strip() for x in l.split('\t')]
            res = reres.match(l[1])
            resnum = int(non_digit.sub('', res.groups()[1]))
            mot = motre.match(l[0])
            if mot:
                start = resnum - 7 + len(mot.groups()[0])
                end = resnum + 7 - len(mot.groups()[2])
                instance = l[0].replace('-', '').upper()
            else:
                start = None
                end = None
                instance = l[0]
            result.append({
                'instance': instance,
                'kinase': l[2],
                'resaa': res.groups()[0],
                'resnum': resnum,
                'score': float(non_digit.sub('', l[3])),
                'substrate': substrate,
                'start': start,
                'end': end
            })
    return result


def pnetworks_interactions():
    result = []
    data = get_phosphonetworks()
    for l in data:
        result.append((l['kinase'], l['substrate']))
    return [list(x) for x in list(set(result))]


def get_depod(organism = 9606):
    
    result = []
    
    reunip = re.compile(r'uniprotkb:([A-Z0-9]+)')
    reptm = re.compile(r'([A-Z][a-z]{2})-([0-9]+)')
    repmidsep = re.compile(r'[,|]\s?')
    
    url = urls.urls['depod']['urls'][0]
    c = curl.Curl(url, silent = False, encoding = 'ascii')
    data = c.result
    data = [x.split('\t') for x in data.split('\n')]
    del data[0]
    
    url_mitab = urls.urls['depod']['urls'][1]
    c_mitab = curl.Curl(url_mitab, silent = False, encoding = 'iso-8859-1')
    data_mitab = c_mitab.result
    data_mitab = [x.split('\t') for x in data_mitab.split('\n')]
    del data_mitab[0]
    
    for i, l in enumerate(data):
        
        if (
            len(l) > 6 and
            l[2] == 'protein substrate' and
            taxonomy.ensure_ncbi_tax_id(
                l[3].split('(')[0].strip()
            ) == organism and
            l[4].strip() != 'N/A'
        ):
            
            enzyme_uniprot = reunip.search(data_mitab[i][0]).groups()[0]
            substrate_uniprot = reunip.search(data_mitab[i][1]).groups()[0]
            
            for enzyme_up, substrate_up in itertools.product(
                    mapping.map_name(
                        enzyme_uniprot,
                        'uniprot',
                        'uniprot'
                    ),
                    mapping.map_name(
                        substrate_uniprot,
                        'uniprot',
                        'uniprot'
                    ),
                ):
                
                for resaa, resnum in reptm.findall(l[4]):
                    
                    resnum = int(resnum)
                    resaa = (
                        common.aminoa_3_to_1_letter[resaa]
                            if resaa in common.aminoa_3_to_1_letter else
                        resaa
                    )
                    
                    result.append({
                        'instance': None,
                        'kinase': enzyme_up,
                        'resaa': resaa,
                        'resnum': resnum,
                        'references': repmidsep.split(l[6].strip()),
                        'substrate': substrate_up,
                        'start': None,
                        'end': None,
                        'typ': 'dephosphorylation',
                    })
    
    return result


def phosphopoint_interactions():
    
    interactions = []
    
    url = urls.urls['phosphopoint']['url']
    
    c = curl.Curl(url, silent = False, large = True)
    
    _ = next(c.result)
    
    return [
        l.strip().split(';')
        for l in c.result
    ]


def phosphopoint_directions():
    
    return [
        l[:2] for l in phosphopoint_interactions()
    ]


def get_acsn():
    
    greek = {
        '_alpha_': 'A',
        '_beta_': 'B',
        '_gamma_': 'C',
        '_delta_': 'D',
        '_epsilon_': 'E'
    }
    regreek = re.compile(r'\b(' + '|'.join(greek.keys()) + r')\b')
    result = []
    url = urls.urls['acsn']['sif']
    c = curl.Curl(url, silent = False)
    data = c.result
    data = [
        x.split('\t')
        for x in data.replace('\r', '').replace('*', '').strip().split('\n')
    ]
    
    for l in data:
        l[0] = regreek.sub('', l[0]).split('_')[0].split('~')[0]
        l[2] = regreek.sub('', l[2]).split('_')[0].split('~')[0]
    
    return data


def get_abs():

    result = []
    url = urls.urls['abs']['url']
    c = curl.Curl(url, silent = False)
    data = c.result
    data = [[x.replace('*', '') for x in xx.split('\t')]
            for xx in data.split('\n')]
    for d in data:
        if len(d) > 2:
            result.append([d[2], d[0]])
    return result


def get_pazar():

    url = urls.urls['pazar']['url_rescued']
    c = curl.Curl(url, silent = False)
    data = c.result
    return [
        list(map(x.split('\t').__getitem__, (1, 4, 10)))
        for x in ''.join(data.values()).split('\n') if len(x) > 0
    ]


def get_htri():
    
    HTRIInteraction = collections.namedtuple(
        'HTRIInteraction',
        [
            'entrez_tf',
            'genesymbol_tf',
            'entrez_target',
            'genesymbol_target',
            'pubmed',
        ]
    )
    
    c = curl.Curl(
        urls.urls['htri']['url'],
        init_url = urls.urls['htri']['init_url'],
        silent = False,
        follow = False,
        large = True,
    )
    
    data = c.result
    _ = next(c.result)
    
    return [
        HTRIInteraction(
            entrez_tf = fields[1],
            genesymbol_tf = fields[2],
            entrez_target = fields[3],
            genesymbol_target = fields[4],
            pubmed = fields[6],
        )
        for fields in
        (
            row.split(';') for row in data if row.strip()
        )
    ]


def get_oreganno_old(organism = 9606):

    taxids = common.swap_dict(taxonomy.taxids)

    if organism in taxids:
        organism = taxids[organism]

    nsep = re.compile(r'([-A-Za-z0-9]{3,})[\s/\(]*.*')
    nrem = re.compile(r'[-/]')
    result = []
    url = urls.urls['oreganno_old']['url']
    c = curl.Curl(url, silent = False)
    data = c.result
    data = [[xx.strip() for xx in x.split('\t')] for x in data.split('\n')
            if len(x) > 0][1:]
    for l in data:
        if l[0] == organism and \
            l[10] == 'TRANSCRIPTION FACTOR BINDING SITE' and \
            l[3] == 'POSITIVE OUTCOME' and \
                not l[11].startswith('UNKNOWN') and not l[14].startswith('UNKNOWN'):
            result.append([
                l[14]
                if len(l[14]) < 3 else nrem.sub('',
                                                nsep.findall(l[14])[0]), l[11]
                if len(l[11]) < 3 else nrem.sub('', nsep.findall(l[11])[0]),
                l[18]
            ])
    return result

def get_oreganno(organism = 9606):

    taxids = taxonomy.phosphoelm_taxids

    if organism in taxids:
        organism = taxids[organism]

    nsep = re.compile(r'([-A-Za-z0-9]{3,})[\s/\(]*.*')
    nrem = re.compile(r'[-/]')
    result = []

    url = urls.urls['oreganno']['url']
    c = curl.Curl(url, silent = False, large = True)
    data = c.result
    _ = next(data)

    for l in data:

        if not l:

            continue

        l = [x.strip() for x in l.split('\t')]

        if (l[1] == organism and
            l[3] == 'TRANSCRIPTION FACTOR BINDING SITE' and
            l[2] == 'POSITIVE OUTCOME' and
            l[4] != 'N/A' and
            l[7] != 'N/A'):

            yield (
                l[7]
                if len(l[7]) < 3 else nrem.sub('',
                                                nsep.findall(l[7])[0]), l[4]
                if len(l[4]) < 3 else nrem.sub('', nsep.findall(l[4])[0]),
                l[11] if l[11] != 'N/A' else ''
            )

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


def get_pathwaycommons(sources = None, types = None, sources_separated = True):

    result = {}
    interactions = []

    if type(types) is list:
        types = set(types)

    source_names = {
        'wp': 'WikiPathways',
        'kegg': 'KEGG',
        'bind': 'BIND',
        'intact': 'IntAct',
        'intact_complex': 'IntAct',
        'panther': 'PANTHER',
        'pid': 'NCI-PID',
        'reactome': 'Reactome',
        'dip': 'DIP',
        'hprd': 'HPRD',
        'inoh': 'INOH',
        'netpath': 'NetPath',
        'biogrid': 'BioGRID',
        'corum': 'CORUM',
        'psp': 'PhosphoSite'
    }

    directed = set([
        'state-change', 'controls-state-change-of', 'controls-transport-of',
        'controls-phosphorylation-of'
    ])

    sources = list(source_names.keys()) \
        if sources is None else sources

    prg = progress.Progress(
        len(sources), 'Processing PathwayCommons', 1, percent = False)

    url = urls.urls['pwcommons']['url']

    for s in sources:

        prg.step()
        surl = url % s
        c = curl.Curl(surl, silent = False, large = True)

        for l in c.result:

            if hasattr(l, 'decode'):

                l = l.decode('ascii')

            l = l.strip().split('\t')

            if types is None or l[1] in types:

                if sources_separated:
                    l.append(source_names[s])
                    interactions.append(l)

                else:
                    pair = (l[0], l[2])

                    if pair not in result:

                        result[pair] = [set([]), set([]), 0]

                    result[pair][0].add(source_names[s])
                    result[pair][1].add(l[1])

                    if l[1] in directed:
                        result[pair][2] = 1

    if not sources_separated:

        for pair, details in iteritems(result):

            interactions.append([
                pair[0], pair[1], ';'.join(details[0]), ';'.join(details[1]),
                str(details[2])
            ])

    return interactions


def go_annotations_uniprot(organism = 9606, swissprot = 'yes'):
    """
    Deprecated, should be removed soon.
    """

    rev = '' if swissprot is None \
        else ' AND reviewed:%s' % swissprot
    query = 'organism:%u%s' % (int(organism), rev)
    url = urls.urls['uniprot_basic']['url']
    get = {'query': query, 'format': 'tab', 'columns': 'id,go-id'}
    c = curl.Curl(url, get = get, silent = False)
    data = c.result
    return dict([(x[0], [go.strip() for go in x[1].split(';')])
                 for x in [x.split('\t') for x in data.split('\n')]
                 if len(x) > 1])


def go_annotations_goa(organism = 'human'):
    """
    Downloads GO annotation from UniProt GOA.
    """

    organism = (
        taxonomy.taxids[organism]
            if isinstance(organism, int) else
        organism
    )

    annot = dict(
        (asp, collections.defaultdict(set))
        for asp in ('C', 'P', 'F')
    )

    url = urls.urls['goa']['ebi_url'] % (organism.upper(), organism)
    c = curl.Curl(url, silent = False, large = True)

    for line in c.result:

        if not line or line[0] == '!':
            continue

        line = line.strip().split('\t')

        annot[line[8]][line[1]].add(line[4])

    return annot


# synonym for the default method
go_annotations = go_annotations_goa


def go_ancestors_goose(aspects = ('C','F','P')):
    """
    Queries the ancestors of GO terms by AmiGO goose.

    Returns dict of sets where keys are GO accessions and values are sets
    of their ancestors.

    :param tuple aspects:
        GO aspects: `C`, `F` and `P` for cellular_component,
        molecular_function and biological_process, respectively.
    """

    aspects_part = ''
    respaces = re.compile(r'[\s\n]+')

    ontologies = {
        'C': 'cellular_component',
        'F': 'molecular_function',
        'P': 'biological_process',
    }

    if set(aspects) != {'C', 'F', 'P'}:

        aspects_part = 'WHERE (%s)' % (
            ' OR '.join(
                'term.term_type = "%s"' % ontologies[asp]
                for asp in aspects
            )
        )

    sql_path = os.path.join(common.DATA, 'goose_ancestors.sql')

    with open(sql_path, 'r') as fp:

        query = fp.read()

    query = query % aspects_part
    query = respaces.sub(r' ', query).strip()

    url = urls.urls['goose']['url'] % query

    c = curl.Curl(url, silent = False, large = True)

    ancestors = collections.defaultdict(set)

    for l in c.result:

        l = l.strip().split('\t')

        ancestors[l[0]].add(l[1])

    return ancestors


def go_ancestors_quickgo(aspects = ('C', 'F', 'P')):
    """
    Queries the ancestors of GO terms by QuickGO REST API.

    Returns dict of sets where keys are GO accessions and values are sets
    of their ancestors.

    :param tuple aspects:
        GO aspects: `C`, `F` and `P` for cellular_component,
        molecular_function and biological_process, respectively.
    """

    desc = go_descendants_quickgo(aspects = aspects)

    return go_descendants_to_ancestors(desc)


# synonym for the default method
go_ancestors = go_ancestors_quickgo


def go_descendants_to_ancestors(desc):
    """
    Turns a dict of descendants to dict of ancestors by swapping the
    relationships. This way descendants will be the keys and their ancestors
    will be the values.
    """

    ancestors = {}

    for asp, dct in iteritems(desc):

        ancestors[asp] = collections.defaultdict(set)

        for anc_term, des in iteritems(dct):

            for des_term, rel in des:

                ancestors[asp][des_term].add((anc_term, rel))

    return ancestors


def go_descendants_goose(aspects = ('C','F','P')):
    """
    Queries descendants of GO terms by AmiGO goose.

    IMPORTANT:
    This is not the preferred method any more to get descendants.
    Recently the preferred method to access GO annotations is
    ``pypath.dataio.go_descendants_quickgo()``.
    The data in GO MySQL instances has not been updated since Dec 2016.
    Unfortunately the providers ceased to support MySQL, the most flexible
    and highest performance access to GO data. The replacement is Solr
    which is far from providing the same features as MySQL, for example
    it is unable to provide GO graph relationships. Other service is QuickGO
    which is up to date and has nice ways to query the ontology.

    Returns dict of sets where keys are GO accessions and values are sets
    of their descendants.

    :param tuple aspects:
        GO aspects: `C`, `F` and `P` for cellular_component,
        molecular_function and biological_process, respectively.
    """

    desc = collections.defaultdict(set)

    anc = go_ancestors_goose(aspects = aspects)

    for term, ancs in iteritems(anc):

        for terma in ancs:

            desc[terma].add(term)

    return desc


def go_descendants_quickgo(
        aspects = ('C', 'F', 'P'),
        terms = None,
        relations = None,
        quickgo_download_size = 500,
    ):
    """
    Queries descendants of GO terms by QuickGO REST API.

    Returns dict of sets where keys are GO accessions and values are sets
    of their descendants.

    :param tuple aspects:
        GO aspects: `C`, `F` and `P` for cellular_component,
        molecular_function and biological_process, respectively.
    :param dict terms:
        Result from ``go_terms_solr``. If ``None`` the method will be called.
    """


    def download_in_chunks(terms, chunk_size, target = None):

        target = target or collections.defaultdict(set)

        paginator = common.paginate(terms, chunk_size)

        for p, terms_part in enumerate(paginator):

            url = urls.urls['quickgo_rest']['desc'] % (
                ','.join(terms_part),
                '?relations = %s' % relations_part,
            )

            c = curl.Curl(
                url,
                req_headers = req_headers,
                silent = True,
                large = True,
            )

            try:
                result = json.load(c.fileobj)
            except json.decoder.JSONDecodeError:
                done = chunk_size * p
                remaining = terms[done:]
                new_chunk_size = chunk_size // 2

                if new_chunk_size < 10:

                    _log(
                        'Failed to download QuickGO, tried to decrease the '
                        'number of terms in each query, went below 10 terms '
                        'per query but still getting erroneous JSON. '
                        'This might be due to very slow network connection. '
                        'You might increase the timeout of CURL. '
                        'But then it will take forever.'
                    )

                    return target

                return download_in_chunks(
                    terms = remaining,
                    chunk_size = new_chunk_size,
                    target = taret,
                )

            for res in result['results']:

                if 'children' not in res:

                    continue

                target[res['id']].update(
                    set(
                        (child['id'], child['relation'])
                        for child in res['children']
                    )
                )

        return target


    desc = {}

    terms = terms or go_terms_quickgo(aspects = aspects)
    relations = relations or ('is_a', 'part_of', 'occurs_in', 'regulates',)

    req_headers = ['Accept:application/json']

    relations_part = ','.join(relations)

    for asp in aspects:

        desc[asp] = download_in_chunks(
            terms = list(terms[asp].keys()),
            chunk_size = quickgo_download_size,
        )

    return desc


# synonym for the default method
go_descendants = go_descendants_quickgo


def go_terms_solr(aspects = ('C', 'F', 'P')):
    """
    Queries GO terms by AmiGO Solr.

    Returns dict of dicts where upper level keys are one letter codes of the
    aspects `C`, `F` and `P` for cellular_component, molecular_function and
    biological_process, respectively. Lower level keys are GO accessions
    and values are names of the terms.

    :param tuple aspects:
        GO aspects: `C`, `F` and `P` for cellular_component,
        molecular_function and biological_process, respectively.
    """

    reamp = re.compile(r'[\s\n\r]+([&\?])')
    relin = re.compile(r'[\s\n\r]+')

    ontologies = {
        'C': 'cellular_component',
        'F': 'molecular_function',
        'P': 'biological_process',
    }
    ontol_short = dict(reversed(i) for i in ontologies.items())

    terms = dict((a, {}) for a in aspects)

    query = '''
        ?q = document_category:"ontology_class" AND
            idspace:GO AND
            is_obsolete:0
        &rows = 9999999
        &start = 0
        &fl = annotation_class,annotation_class_label,source
    '''

    query = relin.sub(' ', reamp.sub(r'\1', query.strip()))

    # downloading data
    url = urls.urls['golr']['url'] % query

    c = curl.Curl(url, silent = False, large = True)

    # parsing XML by lxml.etree.iterparse
    parser = etree.iterparse(c.fileobj, events = ('start', 'end'))
    root = next(parser)
    used_elements = []

    for ev, elem in parser:

        if ev == 'end' and elem.tag == 'doc':

            asp  = elem.find('.//str[@name="source"]').text
            asp  = ontol_short[asp]

            if asp not in aspects:

                continue

            term = elem.find('.//str[@name="annotation_class"]').text
            name = elem.find('.//str[@name="annotation_class_label"]').text

            terms[asp][term] = name

        used_elements.append(elem)

        # removing used elements to keep memory low
        if len(used_elements) > 1000:

            for _ in xrange(500):

                e = used_elements.pop(0)
                e.clear()

    # closing the XML
    c.fileobj.close()
    del c

    return terms


def go_terms_quickgo(aspects = ('C','F','P')):
    """
    Queries GO terms by the QuickGO REST API.

    Return dict of dicts where upper level keys are one letter codes of the
    aspects `C`, `F` and `P` for cellular_component, molecular_function and
    biological_process, respectively. Lower level keys are GO accessions
    and values are names of the terms.

    :param tuple aspects:
        GO aspects: `C`, `F` and `P` for cellular_component,
        molecular_function and biological_process, respectively.
    """

    ontologies = {
        'C': 'cellular_component',
        'F': 'molecular_function',
        'P': 'biological_process',
    }
    ontol_short = dict(reversed(i) for i in ontologies.items())

    result = dict((a, {}) for a in aspects)
    url = urls.urls['quickgo_rest']['terms']
    last_page = 9999999
    this_page = 1
    prg = progress.Progress(
        name = 'Downloading data from QuickGO',
        interval = 1,
    )

    while this_page <= last_page:

        page_url = url % this_page

        c = curl.Curl(page_url, silent = True)

        this_result = json.loads(c.result)
        last_page = this_result['pageInfo']['total']


        for res in this_result['results']:

            if 'aspect' not in res:

                continue

            asp = ontol_short[res['aspect']]

            if res['isObsolete'] or asp not in aspects:

                continue

            result[asp][res['id']] = res['name']

        if prg.total is None:

            prg.set_total(last_page)

        prg.step()

        this_page += 1

    return result


# synonym for the default method
go_terms = go_terms_quickgo


def go_terms_goose(aspects = ('C','F','P')):
    """
    Queries GO terms by AmiGO goose.

    Return dict of dicts where upper level keys are one letter codes of the
    aspects `C`, `F` and `P` for cellular_component, molecular_function and
    biological_process, respectively. Lower level keys are GO accessions
    and values are names of the terms.

    :param tuple aspects:
        GO aspects: `C`, `F` and `P` for cellular_component,
        molecular_function and biological_process, respectively.
    """

    aspects_part = ''
    respaces = re.compile(r'[\s\n]+')

    ontologies = {
        'C': 'cellular_component',
        'F': 'molecular_function',
        'P': 'biological_process',
    }
    ontol_short = dict(reversed(i) for i in ontologies.items())

    if set(aspects) != {'C', 'F', 'P'}:

        aspects_part = 'WHERE (%s)' % (
            ' OR '.join(
                'term.term_type = "%s"' % ontologies[asp]
                for asp in aspects
            )
        )

    sql_path = os.path.join(common.DATA, 'goose_terms.sql')

    with open(sql_path, 'r') as fp:

        query = fp.read()

    query = query % aspects_part
    query = respaces.sub(r' ', query).strip()

    url = urls.urls['goose']['url'] % query

    c = curl.Curl(url, silent = False, large = True)

    terms = {'P': {}, 'C': {}, 'F': {}}

    for l in c.result:

        l = l.strip().split('\t')

        if l[1] not in ontol_short:

            continue

        aspect = ontol_short[l[1]]

        terms[aspect][l[2]] = l[0]

    return terms


def go_annotations_quickgo(
        organism = 9606,
        aspects = ('C','F','P'),
        relations = ('is_a', 'part_of'),
    ):
    """
    Queries GO annotations by QuickGO REST API.

    IMPORTANT:
    Recently the preferred method to access GO annotations is
    ``pypath.dataio.go_annotations_goa()``.
    Contrary to its name QuickGO is super slow, otherwise it should yield
    up to date data, identical to the GOA file.

    Returns terms in dict of dicts and annotations in dict of dicts of sets.
    In both dicts the keys are aspects by their one letter codes.
    In the term dicts keys are GO accessions and values are their names.
    In the annotation dicts keys are UniProt IDs and values are sets
    of GO accessions.

    :param int organism:
        NCBI Taxonomy ID of one organism. Default is human (9606).
    :param tuple aspects:
        GO aspects: `C`, `F` and `P` for cellular_component,
        molecular_function and biological_process, respectively.
    :param list uniprots:
        Optionally a list of UniProt IDs. If `None`, results for all proteins
        returned.
    """

    annot = dict((a, collections.defaultdict(set)) for a in aspects)

    ontologies = {
        'C': 'cellular_component',
        'F': 'molecular_function',
        'P': 'biological_process',
    }
    ontol_short = dict(reversed(i) for i in ontologies.items())

    url = urls.urls['quickgo_rest']['annot']

    aspects_part = ','.join(ontologies[a] for a in aspects)
    relations_part = ','.join(relations)

    req_headers = ['Accept:text/tsv']

    page = 1

    while True:

        this_url = url % (
            aspects_part, # aspect
            relations_part, # goUsageRelationships
            organism, # taxonId
            page,
        )

        c = curl.Curl(
            url = this_url,
            req_headers = req_headers,
            silent = False,
            large = True
        )

        _ = next(c.result) # the header row

        for line in c.result:

            line = line.strip().split('\t')

            if line[3] not in relations:

                continue

            annot[line[5]][line[1]].add(line[4])

        page += 1

    return annot


def go_annotations_solr(
        organism = 9606,
        aspects = ('C', 'F', 'P'),
        references = False,
    ):
    """
    Queries GO annotations by AmiGO Solr.

    Before other methods have been provided to access GO.
    Now this is the preferred method to get annotations.
    Returns terms in dict of dicts and annotations in dict of dicts of sets.
    In both dicts the keys are aspects by their one letter codes.
    In the term dicts keys are GO accessions and values are their names.
    In the annotation dicts keys are UniProt IDs and values are sets
    of GO accessions.

    :param int organism:
        NCBI Taxonomy ID of one organism. Default is human (9606).
    :param tuple aspects:
        GO aspects: `C`, `F` and `P` for cellular_component,
        molecular_function and biological_process, respectively.
    :param bool references:
        Retrieve the references (PubMed IDs) for the annotations.
        Currently not implemented.
    """

    reamp = re.compile(r'[\s\n\r]+([&\?])')
    relin = re.compile(r'[\s\n\r]+')

    annot = dict((a, collections.defaultdict(set)) for a in aspects)

    ontologies = {
        'C': 'cellular_component',
        'F': 'molecular_function',
        'P': 'biological_process',
    }
    ontol_short = dict(reversed(i) for i in ontologies.items())

    # assembling the query

    if len(aspects) < 3:

        aspects_part = ' AND (%s)' % (
            ' OR '.join('aspect:%s' % a for a in aspects)
        )

    else:

        aspects_part = ''

    refs_part = ',reference' if references else ''

    query = '''
        ?q = taxon:"NCBITaxon:%u" AND
            type:protein AND
            document_category:annotation AND
            source:UniProtKB%s
        &rows = 9999999
        &start = 0
        &fl = bioentity,annotation_class,aspect%s
    ''' % (
        organism,
        aspects_part,
        refs_part
    )

    query = relin.sub(' ', reamp.sub(r'\1', query.strip()))

    # downloading data
    url = urls.urls['golr']['url'] % query
    c = curl.Curl(url, silent = False, large = True)

    # parsing XML by lxml.etree.iterparse
    parser = etree.iterparse(c.fileobj, events = ('start', 'end'))
    root = next(parser)
    used_elements = []

    for ev, elem in parser:

        if ev == 'end' and elem.tag == 'doc':

            id_ = elem.find('.//str[@name="bioentity"]').text

            if not id_.startswith('UniProtKB:'):

                continue

            asp  = elem.find('.//str[@name="aspect"]').text

            if asp not in aspects:

                continue

            term = elem.find('.//str[@name="annotation_class"]').text
            id_  = id_[10:] # removing the `UniProtKB:` prefix

            # adding the term to the annotation dict
            annot[asp][id_].add(term)

        used_elements.append(elem)

        # removing used elements to keep memory low
        if len(used_elements) > 1000:

            for _ in xrange(500):

                e = used_elements.pop(0)
                e.clear()

    # closing the XML
    c.fileobj.close()
    del c

    return terms, annot


def go_annotations_goose(organism = 9606, aspects = ('C','F','P'), uniprots = None):
    """
    Queries GO annotations by AmiGO goose.

    IMPORTANT:
    This is not the preferred method any more to get terms and annotations.
    Recently the preferred method to access GO annotations is
    ``pypath.dataio.go_annotations_solr()``.
    The data in GO MySQL instances has not been updated since Dec 2016.
    Unfortunately the providers ceased to support MySQL, the most flexible
    and highest performance access to GO data. The replacement is Solr
    which is far from providing the same features as MySQL.

    Returns terms in dict of dicts and annotations in dict of dicts of sets.
    In both dicts the keys are aspects by their one letter codes.
    In the term dicts keys are GO accessions and values are their names.
    In the annotation dicts keys are UniProt IDs and values are sets
    of GO accessions.

    :param int organism:
        NCBI Taxonomy ID of one organism. Default is human (9606).
    :param tuple aspects:
        GO aspects: `C`, `F` and `P` for cellular_component,
        molecular_function and biological_process, respectively.
    :param list uniprots:
        Optionally a list of UniProt IDs. If `None`, results for all proteins
        returned.
    """

    aspects_part = ''
    uniprot_part = ''
    respaces = re.compile(r'[\s\n]+')

    ontologies = {
        'C': 'cellular_component',
        'F': 'molecular_function',
        'P': 'biological_process',
    }
    ontol_short = dict(reversed(i) for i in ontologies.items())

    if set(aspects) != {'C', 'F', 'P'}:

        aspects_part = '(%s) AND' % (
            ' OR '.join(
                'term.term_type="%s"' % ontologies[asp]
                for asp in aspects
            )
        )

    if uniprots is not None:

        uniprot_part = 'dbxref.xref_key IN (%s) AND' % (
            ','.join('"%s"' % uniprot for uniprot in uniprots)
        )

    sql_path = os.path.join(common.DATA, 'goose_annotations.sql')

    with open(sql_path, 'r') as fp:

        query = fp.read()

    query = query % (organism, aspects_part, uniprot_part)
    query = respaces.sub(r' ', query).strip()

    url = urls.urls['goose']['url'] % query

    c = curl.Curl(url, silent = False, large = True)

    terms = {'P': {}, 'C': {}, 'F': {}}
    annot = {
        'C': collections.defaultdict(set),
        'F': collections.defaultdict(set),
        'P': collections.defaultdict(set),
    }

    for l in c.result:

        l = l.strip().split('\t')

        aspect = ontol_short[l[1]]

        terms[aspect][l[2]] = l[0]
        annot[aspect][l[5]].add(l[2])

    return terms, annot


def get_go_desc(go_ids, organism = 9606):
    """
    Deprecated, should be removed soon.
    """

    go_ids = (
        ','.join(sorted(go_ids))
        if type(go_ids) in {list, tuple, set} else
        go_ids
    )

    url = urls.urls['quickgo_desc']['url'] % (organism, go_ids)

    c = curl.Curl(
        url, silent = False, large = True, req_headers = {'Accept': 'text/tsv'}
    )
    _ = c.result.readline()

    return set(l.split('\t')[1] for l in c.result)


def get_go_quick(
        organism = 9606,
        slim = False,
        names_only = False,
        aspects = ('C', 'F', 'P'),
    ):
    """
    Deprecated, should be removed soon.

    Loads GO terms and annotations from QuickGO.
    Returns 2 dicts: `names` are GO terms by their IDs,
    `terms` are proteins GO IDs by UniProt IDs.
    """

    ontologies = {
        'C': 'cellular_component',
        'F': 'molecular_function',
        'P': 'biological_process',
    }

    terms = {
        'C': collections.defaultdict(set),
        'F': collections.defaultdict(set),
        'P': collections.defaultdict(set),
    }
    names = {}
    aspects_param = ','.join(sorted(ontologies[a] for a in aspects))
    url = urls.urls['quickgo']['url'] % (
        organism,
        aspects_param,
        '&goUsage = slim' if slim else '',
    )

    c = curl.Curl(url, silent = False, large = True)
    _ = next(c.result)

    for l in result:

        l = l.split('\t')

        if not names_only:

            terms[l[5]][l[1]].add(l[4])

    return {'terms': terms, 'names': names}


def get_goslim(url = None):
    rego = re.compile(r'GO:[0-9]{7}')
    url = url if type(url) in [str, unicode] \
        else urls.urls['goslim_gen']['url']
    c = curl.Curl(url, silent = False)
    data = c.result
    result = []
    for l in data.split('\n'):
        if l.startswith('id:'):
            result += rego.findall(l)
    return result


def netpath_names():

    repwnum = re.compile(r'_([0-9]+)$')
    result = {}
    url = urls.urls['netpath_names']['url']
    c = curl.Curl(url, silent = False)
    html = c.result
    soup = bs4.BeautifulSoup(html, 'html.parser')
    for a in soup.find_all('a'):
        if a.attrs['href'].startswith('pathways'):
            num = repwnum.findall(a.attrs['href'])[0]
            name = a.text
            result[num] = name.strip()
    return result


def netpath_pathway_annotations():

    NetpathPathway = collections.namedtuple(
        'NetpathPathway', ['pathway'],
    )


    result = collections.defaultdict(set)

    url_template = urls.urls['netpath_pw']['url']

    url_main = urls.urls['netpath_pw']['mainpage']
    c = curl.Curl(url_main, cache = False)
    cookie = [
        h.decode().split(':')[1].split(';')[0].strip()
        for h in c.resp_headers
        if h.startswith(b'Set-Cookie')
    ]
    cookie_hdr = ['Cookie: %s' % '; '.join(cookie)]

    pathway_ids = netpath_names()

    for _id, pathway in iteritems(pathway_ids):

        url = url_template % int(_id)
        c = curl.Curl(
            url,
            req_headers = cookie_hdr,
            silent = False,
            encoding = 'iso-8859-1',
        )

        soup = bs4.BeautifulSoup(c.result, 'html.parser')

        for tbl in soup.find_all('table'):

            hdr = tbl.find('td', {'class': 'barhead'})

            if not hdr or not hdr.text.strip().startswith('Molecules Invol'):

                continue

            for td in tbl.find_all('td'):

                genesymbol = td.text.strip()

                if not genesymbol:

                    continue

                uniprots = mapping.map_name(
                    genesymbol,
                    'genesymbol',
                    'uniprot',
                )
                
                for uniprot in uniprots:
                    
                    result[uniprot].add(
                        NetpathPathway(
                            pathway = pathway
                        )
                    )

    return result


def netpath_interactions():
    
    result = []
    repwnum = re.compile(r'NetPath_([0-9]+)_')
    mi = '{net:sf:psidev:mi}'
    url = urls.urls['netpath_psimi']['url']
    c = curl.Curl(url, silent = False)
    
    data = c.result
    data = dict([(k, v) for k, v in iteritems(data) if k.endswith('xml')])
    pwnames = netpath_names()
    
    for pwfile, rawxml in iteritems(data):
        try:
            pwnum = repwnum.findall(pwfile)[0]
        except:
            sys.stdout.write('Error at processing file:\n')
            sys.stdout.write(pwfile)
            sys.stdout.write('\n')
            sys.stdout.flush()
        pwname = pwnames[pwnum]
        root = ET.fromstring(rawxml)
        for e in root.findall(mi + 'entry'):
            thisInt = ()
            db = [
                pr.find(mi + 'primaryRef').attrib['db']
                for pr in e.find(mi + 'source').findall(mi + 'xref')
            ]
            refs = []
            mets = []
            for ex in e.find(mi + 'experimentList').findall(
                    mi + 'experimentDescription'):
                for pm in ex.find(mi + 'bibref').iter(mi + 'primaryRef'):
                    if pm.attrib['db'] == 'pubmed':
                        refs.append(pm.attrib['id'])
                for me in ex.find(mi + 'interactionDetectionMethod').\
                        iter(mi + 'shortLabel'):
                    mets.append(me.text)
            mols = {}
            for mo in e.find(mi + 'interactorList').findall(mi + 'interactor'):
                iid = mo.attrib['id']
                name = mo.find(mi + 'names').find(mi + 'shortLabel').text
                entrez = ''
                if mo.find(mi + 'xref') is not None:
                    entrez = ';'.join([
                        ac.attrib['id']
                        for ac in mo.find(mi + 'xref')
                        .findall(mi + 'secondaryRef')
                        if ac.attrib['db'] == 'Entrez gene'
                    ])
                mols[iid] = (name, entrez)
            theInt = e.find(mi + 'interactionList').find(mi + 'interaction')
            for p in theInt.find(mi + 'participantList').findall(
                    mi + 'participant'):
                pid = p.find(mi + 'interactorRef').text
                roles = ''
                if p.find(mi + 'experimentalRoleList') is not None:
                    roles = ';'.join([
                        rl.find(mi + 'names').find(mi + 'shortLabel').text
                        for rl in p.find(mi + 'experimentalRoleList')
                        .findall(mi + 'experimentalRole')
                    ])
                mols[pid] += (roles, )
            intTyp = theInt.find(mi + 'interactionType').find(mi + 'names')\
                .find(mi + 'shortLabel').text
            molkeys = list(mols.keys())
            for i in range(0, len(mols) - 1):
                for j in range(i, len(mols)):
                    A = mols[molkeys[i]][0:2]
                    B = mols[molkeys[j]][0:2]
                    result.append(
                        list(A) + list(B) +
                        [';'.join(refs), ';'.join(mets), intTyp, pwname])
    return result


def get_pubmeds(pmids):
    pmids = [str(pmid) for pmid in pmids]
    url = urls.urls['pubmed-eutils']['url']
    cache = len(pmids) < 10
    data = {}
    prg = progress.Progress(
        len(pmids) / 100 + 1,
        'Retrieving data from NCBI e-utils',
        1,
        percent = False)
    for offset in xrange(0, len(pmids), 100):
        prg.step()
        post = {
            'id': ','.join(pmids[offset:offset + 100]),
            'retmode': 'json',
            'db': 'pubmed'
        }
        for i in xrange(3):
            try:
                c = curl.Curl(
                    url,
                    silent = False,
                    cache = cache,
                    post = post,
                    override_post = True,
                )
                res = c.result
                data = dict([(k, v)
                             for k, v in iteritems(json.loads(res)['result'])]
                            + [(k, v) for k, v in iteritems(data)])
                break
            except ValueError:
                sys.stdout.write('\t:: Error in JSON, retry %u\n' % i)
                sys.stdout.flush()
    prg.terminate()
    return data


def get_lincs_compounds():
    sys.stdout.write(
        '\n\tReturned dict has names, brand names or company specific\n'
        '\tIDs of compounds as keys, and tuples of PubChem, ChEMBL, ChEBI, InChi, \n'
        '\tInChi Key, SMILES and LINCS as values.\n\n')
    sys.stdout.flush()
    c = curl.Curl(urls.urls['lincs-compounds']['url'], silent = False)
    return dict(
        [(key, pair[1])
         for pair in [([
             it for sl in [
                 filter(lambda z: len(z) > 0, y.split(';')) for y in x[1:4]
                 if len(y) > 0
             ] for it in sl
         ], (x[4], '' if len(x[7]) == 0 else 'CHEMBL%s' % x[7], ''
             if len(x[8]) == 0 else 'CHEBI:%s' % x[8], x[9], x[10], x[11], x[3]
             )) for x in [[b.strip() for b in a.split('\t')] for a in ''.join([
                 s.replace(',', '\t') if i % 2 == 0 else s.replace('\n', '')
                 for i, s in enumerate(c.result.split('"'))
             ]).split('\n')[1:] if len(a) > 0]] for key in pair[0]])


def ramilowski_interactions(putative = False):
    """
    Downloads and processes ligand-receptor interactions from
    Supplementary Table 2 of Ramilowski 2015.
    """

    c = curl.Curl(urls.urls['rami']['url'], silent = False, large = True)
    xlsname = c.fname
    del(c)
    raw = inputs_common.read_xls(xlsname, 'All.Pairs')[1:]

    return [
        [
            r[1],
            r[3],
            r[13].replace(' ', ''),
            ';'.join(filter(len, itertools.chain(r[5:11], [r[15]])))
        ]
        for r in raw
        if r[15] != 'EXCLUDED not ligand' and (
            putative or r[15] != 'putative'
        )
    ]

    return raw

def kirouac2010_interactions():

    rename = re.compile(r'[A-Z]{2}[A-Z0-9][-A-Z0-9]*')
    rerange = re.compile(r'([0-9])-([0-9])')
    reslash = re.compile(r'.*?([A-Z0-9]{1,3}/[/A-Z0-9]+)')

    def get_names(s):

        names = set([])
        prev = None

        for n in s.split():

            m = rename.findall(n)

            if m:

                prev = m

                m = reslash.match(n)

                if m:

                    for post in m.groups()[0].split('/'):

                        for pre in prev:

                            names.add('%s%s' % (pre, post))

                else:

                    m = rerange.match(n)

                    if m:

                        intv = m.groups()

                        for post in range(int(intv[0]), int(intv[1]) + 1):

                            for pre in prev:

                                names.add('%s%u' % (pre, post))

                    else:

                        names.update(prev)

            prev = None

        return names
    
    
    init_url = urls.urls['kirouac2010']['init_url']
    req_headers = [
        (
            'User-Agent: Mozilla/5.0 (X11; Linux x86_64; rv:68.0) '
            'Gecko/20100101 Firefox/68.0'
        ),
    ]
    url = urls.urls['kirouac2010']['url']
    
    c00 = curl.Curl(url, call = False, process = False)
    
    if (
        not os.path.exists(c00.cache_file_name) or
        os.path.getsize(c00.cache_file_name) == 0
    ):
        
        _log('Kirouac 2010 download: requesting website cookie.')
        
        c0 = curl.Curl(
            init_url,
            silent = True,
            large = False,
            req_headers = req_headers,
            follow = False,
            cache = False,
        )
        
        cookies = []
        
        if hasattr(c0, 'resp_headers'):
            
            for hdr in c0.resp_headers:
                
                if hdr.lower().startswith(b'set-cookie'):
                    
                    cookie = hdr.split(b':')[1].split(b';')[0].strip()
                    
                    if cookie not in cookies:
                        
                        cookies.append(cookie.decode('ascii'))
            
            cookies = '; '.join(cookies)
            
            req_headers.append('Cookie: %s' % cookies)
            
            _log('Response header: %s' % str(c0.resp_headers))
            _log('Cookies: %s' % str(cookies))
            _log('Request header: %s' % str(req_headers))
        
        os.remove(c00.cache_file_name)
    
    c = curl.Curl(
        url,
        silent = False,
        large = True,
        req_headers = req_headers,
    )
    xlsname = c.fname
    del(c)
    tbl = inputs_common.read_xls(xlsname, sheet = 'S12')

    result = []

    for r in tbl[2:]:

        namesA = get_names(r[0])
        namesB = get_names(r[1])

        result.extend(list(itertools.product(namesA, namesB)))

    return result


def get_hpmr_old():
    """
    Deprecated, should be removed soon.

    Downloads and processes the list of all human receptors from
    human receptor census (HPMR -- Human Plasma Membrane Receptome).
    Returns list of GeneSymbols.
    """

    c = curl.Curl(urls.urls['hpmr']['url'], silent = False)
    html = c.result
    soup = bs4.BeautifulSoup(html, 'html.parser')

    gnames = [
        row[1].text
        for row in (
            tr.find_all('td')
            for tr in soup.find(
                'table', {'class': 'gridtable'}
            ).find_all('tr')
        )
        if len(row) > 1 and not row[1].text.lower().startswith('similar')
    ]

    return common.uniq_list(gnames)


def hpmr_interactions_old():
    """
    Deprecated, should be removed soon.

    Downloads ligand-receptor and receptor-receptor interactions from the
    Human Plasma Membrane Receptome database.
    """

    cachefile = os.path.join(
        settings.get('cachedir'),
        'hpmr_interactions',
    )

    if os.path.exists(cachefile):

        with open(cachefile, 'r') as fp:

            result = [r.split('\t') for r in fp.read().split('\n')[1:]]

        return result

    rerecname = re.compile(r'Receptor ([A-z0-9]+) interacts with:')
    reint = re.compile(r'(Receptor|Ligand) ([A-z0-9]+) -')
    rerefid = re.compile(r'list_uids=([- \.:,0-9A-z]+)')

    result = []
    recpages = []

    c = curl.Curl(urls.urls['hpmri']['browse'])
    soup = bs4.BeautifulSoup(c.result, 'html.parser')

    for rec in soup.find_all('a', {'title': 'Open Receptor Page'}):

        recpages.append(rec.attrs['href'])

    prg = progress.Progress(len(recpages), 'Downloading HPMR data', 1)

    for url in recpages:

        prg.step()

        c = curl.Curl(url)

        if c.result is None:

            #print('No receptor page: %s' % url)
            continue

        soup = bs4.BeautifulSoup(c.result, 'html.parser')
        ints = soup.find('div', {'id': 'GeneInts'})

        if not ints:

            #print('No interactions: %s' % url)
            continue

        recname = rerecname.search(
            ints.find_previous_sibling('span').text
        )
        recname = recname.groups()[0] if recname else 'Unknown'

        if recname == 'Unknown':

            # print('Could not find receptor name: %s' % url)
            continue

        for td in ints.find_all('td'):

            interactors = []

            for span in td.find_all('span', {'class': 'IntRow'}):

                ints = reint.search(span.text)

                if ints:

                    interactors.append(ints.groups())

            references = []

            for ref in td.find_all(
                'a', {'title': 'click to open reference in new window'}):

                references.append(
                    rerefid.search(ref.attrs['href']).groups()[0]
                )

            result.extend([
                [recname, i[0], i[1], ';'.join(references)]
                for i in interactors
            ])

    prg.terminate()

    with open(cachefile, 'w') as fp:

        fp.write('%s\n' % '\t'.join([
            'receptor', 'partner', 'partner_type', 'references'
        ]))

        fp.write('\n'.join('\t'.join(r) for r in result))

    return result


def get_hpmr(use_cache = None):
    """
    Downloads ligand-receptor and receptor-receptor interactions from the
    Human Plasma Membrane Receptome database.
    """


    def get_complex(interactors, typ, recname = None, references = None):
        """
        typ : str
            `Receptor` or `Ligand`.
        """

        components = [i[1] for i in interactors if i[0] == typ]

        if typ == 'Receptor' and recname:

            components.append(recname)

        if len(components) == 1:

            return components[0]

        elif len(components) > 1:

            return components


    cachefile = settings.get('hpmr_preprocessed')
    use_cache = (
        use_cache
            if isinstance(use_cache, bool) else
        settings.get('use_intermediate_cache')
    )

    if os.path.exists(cachefile) and use_cache:

        _log('Reading HPMR data from cache file `%s`.' % cachefile)

        return pickle.load(open(cachefile, 'rb'))

    rerecname = re.compile(r'Receptor ([A-z0-9]+) interacts with:')
    reint = re.compile(r'(Receptor|Ligand) ([A-z0-9]+) -')
    rerefid = re.compile(r'list_uids=([- \.:,0-9A-z]+)')
    refamid = re.compile(r'.*FamId=([0-9\.]+)')

    a_family_title = 'Open Family Page'
    a_receptor_title = 'Open Receptor Page'
    a_titles = {a_family_title, a_receptor_title}

    interactions = []
    complex_interactions = []
    families = {}
    recpages = []

    c = curl.Curl(urls.urls['hpmri']['browse'])
    soup = bs4.BeautifulSoup(c.result, 'html.parser')

    this_family = ('0', None)
    this_subfamily = ('0', None)
    this_subsubfamily = ('0', None)

    for a in soup.find_all('a'):

        a_title = a.attrs['title'] if 'title' in a.attrs else None

        if a_title not in a_titles:

            continue

        if a_title == a_family_title:

            family_id = refamid.match(a.attrs['href']).groups()[0]

            if family_id.startswith(this_subfamily[0]):

                this_subsubfamily = (family_id, a.text)

            elif family_id.startswith(this_family[0]):

                this_subfamily = (family_id, a.text)
                this_subsubfamily = ('0', None)

            else:

                this_family = (family_id, a.text)
                this_subfamily = ('0', None)
                this_subsubfamily = ('0', None)

        elif a_title == a_receptor_title:

            recpages.append((
                a.attrs['href'],
                this_family[1],
                this_subfamily[1],
                this_subsubfamily[1],
            ))

    prg = progress.Progress(len(recpages), 'Downloading HPMR data', 1)

    i_complex = 0

    for url, family, subfamily, subsubfamily in recpages:

        prg.step()

        c = curl.Curl(url)

        if c.result is None:

            #print('No receptor page: %s' % url)
            continue

        soup = bs4.BeautifulSoup(c.result, 'html.parser')
        ints = soup.find('div', {'id': 'GeneInts'})

        if not ints:

            #print('No interactions: %s' % url)
            continue

        recname = rerecname.search(
            ints.find_previous_sibling('span').text
        )
        recname = recname.groups()[0] if recname else 'Unknown'

        if recname == 'Unknown':

            # print('Could not find receptor name: %s' % url)
            continue

        recname_u = mapping.map_name0(recname, 'genesymbol', 'uniprot')

        if not recname_u:

            continue

        families[recname_u] = (
            family,
            subfamily,
            subsubfamily,
        )

        for td in ints.find_all('td'):

            interactors = []

            for span in td.find_all('span', {'class': 'IntRow'}):

                ints = reint.search(span.text)

                if ints:

                    interactors.append(ints.groups())

            references = []

            for ref in td.find_all(
                'a', {'title': 'click to open reference in new window'}):

                references.append(
                    rerefid.search(ref.attrs['href']).groups()[0].strip()
                )

            interactors_u = []

            for role, genesymbol in interactors:

                uniprot = (
                    mapping.map_name0(genesymbol, 'genesymbol', 'uniprot')
                )

                if uniprot:

                    interactors_u.append((role, uniprot))

            interactions.extend([
                [recname_u, i[0], i[1], ';'.join(references)]
                for i in interactors_u
            ])

            rec_complex = get_complex(
                interactors_u,
                'Receptor',
                recname = recname_u,
                references = references,
            )
            lig_complex = get_complex(
                interactors_u,
                'Ligand',
                references = references,
            )

            if (
                isinstance(rec_complex, list) or
                isinstance(lig_complex, list)
            ):

                complex_interactions.append((lig_complex, rec_complex))

    prg.terminate()

    result = {
        'interactions': interactions,
        'families': families,
        'complex_interactions': complex_interactions,
    }

    pickle.dump(result, open(cachefile, 'wb'))

    return result


def hpmr_complexes(use_cache = None):

    hpmr_data = get_hpmr(use_cache = use_cache)

    complexes = {}

    i_complex = 0

    for components in itertools.chain(*hpmr_data['complex_interactions']):

        if isinstance(components, list):

            cplex = intera.Complex(
                components = components,
                sources = 'HPMR',
                ids = 'HPMR-COMPLEX-%u' % i_complex,
            )

            complexes[cplex.__str__()] = cplex

    return complexes


def hpmr_interactions(use_cache = None):

    hpmr_data = get_hpmr(use_cache = use_cache)

    return hpmr_data['interactions']


def hpmr_annotations(use_cache = None):

    annot = collections.defaultdict(set)

    HPMRAnnotation = collections.namedtuple(
        'HPMRAnnotation',
        ('role', 'mainclass', 'subclass', 'subsubclass'),
    )

    hpmr_data = get_hpmr(use_cache = use_cache)

    for i in hpmr_data['interactions']:

        # first partner is always a receptor
        # (because ligand pages simply don't work on HPMR webpage)
        args1 = ('Receptor',) + (
            hpmr_data['families'][i[0]]
                if i[0] in hpmr_data['families'] else
            (None, None, None)
        )
        # the second is either a ligand or another receptor
        args2 = (i[1],) + (
            hpmr_data['families'][i[2]]
                if i[2] in hpmr_data['families'] else
            (None, None, None)
        )

        annot[i[0]].add(HPMRAnnotation(*args1))
        annot[i[2]].add(HPMRAnnotation(*args2))

    for uniprot, classes in iteritems(hpmr_data['families']):

        args = ('Receptor',) + classes

        annot[uniprot].add(HPMRAnnotation(*args))

    return dict(annot)


def adhesome_interactions():

    AdhesomeInteraction = collections.namedtuple(
        'AdhesomeInteraction',
        ['source', 'target', 'effect', 'type', 'pmid'],
    )

    url = urls.urls['adhesome']['interactions']

    c = curl.Curl(url, large = True, silent = False)

    data = csv.DictReader(c.result, delimiter = ',')

    result = []

    for rec in data:

        result.append(
            AdhesomeInteraction(
                source = rec['Source'],
                target = rec['Target'],
                effect = rec['Effect'],
                type   = rec['Type'],
                pmid   = rec['PMID'],
            )
        )

    return result




def get_cpad():

    url = urls.urls['cpad']['url']

    c = curl.Curl(url, silent = False, large = True, encoding = 'iso-8859-1')

    reader = csv.DictReader(c.result, delimiter = '\t')

    return reader


def cpad_annotations(include_unknown_type = False):

    CpadAnnotation = collections.namedtuple(
        'CpadAnnotation',
        [
            'regulator_type',
            'effect_on_pathway',
            'pathway',
            'effect_on_cancer',
            'effect_on_cancer_outcome',
            'cancer',
            'pathway_category',
        ]
    )

    cpad = get_cpad()

    result = collections.defaultdict(set)

    for rec in cpad:

        if rec['Regulator'] == 'NULL':

            continue

        for regulator in rec['Regulator'].split(' and '):

            uniprot = mapping.map_name0(regulator, 'genesymbol', 'uniprot')

            if uniprot:

                regulator_name = uniprot
                regulator_type = 'protein'

            else:

                mirbase = mapping.map_name(
                    'hsa-%s' % regulator,
                    'mir-mat-name',
                    'mirbase',
                )

                if not mirbase:

                    mirbase = mapping.map_name(
                        'hsa-%s' % regulator,
                        'mir-name',
                        'mirbase',
                    )

                if mirbase:

                    regulator_name = mirbase
                    regulator_type = 'mirna'

                else:

                    if include_unknown_type:

                        regulator_name = regulator
                        regulator_type = 'unknown'

                    else:

                        continue

            if isinstance(regulator_name, common.basestring):

                regulator_name = (regulator_name,)

            for regulator_name_0 in regulator_name:

                record = CpadAnnotation(
                    regulator_type = regulator_type,
                    effect_on_pathway = rec['Regulator_Type'],
                    effect_on_cancer = rec['Regulation_Type'],
                    effect_on_cancer_outcome = rec['Outcome_Description'],
                    pathway = rec['Pathway'],
                    pathway_category = rec['Pathway_Category'],
                    cancer = rec['Cancer'],
                )

                result[regulator_name_0].add(record)

    return result


def cpad_pathway_cancer():
    """
    Collects only the pathway-cancer relationships. Returns sets of records
    grouped in dicts by cancer and by pathway.
    """

    CpadPathwayCancer = collections.namedtuple(
        'CpadPathwayCancer',
        [
            'pathway',
            'cancer',
            'pathway_category',
            'effect_on_cancer',
            'effect_on_cancer_outcome',
        ]
    )

    cpad = get_cpad()

    by_cancer = collections.defaultdict(set)
    by_pathway = collections.defaultdict(set)

    for rec in cpad:

        record = CpadPathwayCancer(
            pathway = rec['Pathway'],
            cancer = rec['Cancer'],
            pathway_category = rec['Pathway_Category'],
            effect_on_cancer = rec['Regulation_Type'],
            effect_on_cancer_outcome = rec['Outcome_Description'],
        )

        by_cancer[rec['Cancer']].add(record)
        by_pathway[rec['Pathway']].add(record)

    return by_cancer, by_pathway


def get_integrins():
    """
    Returns a set of the UniProt IDs of the human integrins from
    Table 1 of Takada et al 2007 (10.1186/gb-2007-8-5-215).
    """

    url = urls.urls['integrins']['url']

    req_headers = [
        'Host: www.ncbi.nlm.nih.gov',
        'User-Agent: Mozilla/5.0 (X11; Linux x86_64; rv:67.0) '\
            'Gecko/20100101 Firefox/67.0',
        'Accept: text/html,application/xhtml+xml,'
            'application/xml;q=0.9,*/*;q=0.8',
        'Accept-Language: en-US,en;q=0.5',
        'Connection: keep-alive',
        'Upgrade-Insecure-Requests: 1',
        'Pragma: no-cache',
        'Cache-Control: no-cache',
    ]

    c = curl.Curl(
        url, silent = False, req_headers = req_headers, large = True,
    )
    soup = bs4.BeautifulSoup(c.fileobj.read(), 'lxml')

    integrins = []

    rows = soup.find_all('tr')

    for tr in rows[1:]:

        cells = [td for td in tr.find_all('td')]
        integrins.append(cells[-1].text.split('}')[-1])

    return mapping.map_names(integrins, 'uniprot', 'uniprot')


def get_tfcensus(classes = ('a', 'b', 'other')):
    """
    Downloads and processes the list of all known transcription factors from
    TF census (Vaquerizas 2009). This resource is human only.
    Returns set of UniProt IDs.
    """

    ensg = []
    hgnc = []
    reensg = re.compile(r'ENSG[0-9]{11}')
    url = urls.urls['vaquerizas2009']['url']
    c = curl.Curl(url, silent = False, large = True)
    f = c.result

    for l in f:

        if len(l) > 0 and l.split('\t')[0] in classes:

            ensg += reensg.findall(l)
            h = l.split('\t')[5].strip()
            if len(h) > 0:
                hgnc.append(h)

    return (
        set.union(*(
            mapping.map_name(e, 'ensembl', 'uniprot')
            for e in ensg
        )) |
        set.union(*(
            mapping.map_name(h, 'genesymbol', 'uniprot')
            for h in hgnc
        ))
    )


def get_guide2pharma(
        organism = 'human',
        endogenous = True,
        process_interactions = True,
        process_complexes = True,
    ):
    """
    Downloads and processes Guide to Pharmacology data.
    Returns list of dicts.

    @organism : str
        Name of the organism, e.g. `human`.
    @endogenous : bool
        Whether to include only endogenous ligands interactions.
    """

    get_taxid = taxonomy.taxid_from_common_name

    if isinstance(organism, common.basestring):

        try:

            organism = taxonomy.taxid_from_common_name(organism)

        except KeyError:

            organism = None

    positives = {
        'agonist', 'activator', 'potentiation', 'partial agonist',
        'inverse antagonist', 'full agonist', 'activation',
        'irreversible agonist', 'positive',
    }
    negatives = {
        'inhibitor', 'antagonist', 'inhibition', 'irreversible inhibition',
        'inverse agonist', 'negative', 'weak inhibition',
        'reversible inhibition',
    }


    GuideToPharmacologyInteraction = collections.namedtuple(
        'GuideToPharmacologyInteraction',
        [
            'ligand',
            'ligand_id_type',
            'target',
            'target_id_type',
            'target_is_ligand',
            'ligand_organism',
            'target_organism',
            'effect',
            'ligand_location',
            'target_type',
            'ligand_endogenous',
            'pubmed_ids',
        ]
    )


    def is_positive(term):

        return term.lower().strip() in positives


    def is_negative(term):

        return term.lower().strip() in negatives


    interactions = []
    complexes = {}

    url = urls.urls['gtp']['url']

    c = curl.Curl(url, silent = False, large = True, encoding = 'utf-8')

    data = csv.DictReader(c.result)

    #return data

    if organism is not None:

        data = [
            d for d in data
            if (
                get_taxid(d['target_species']) == organism and
                organism in set(
                    get_taxid(t)
                    for t in d['ligand_species'].split('|')
                )
            )
        ]

    if endogenous:

        data = [d for d in data if d['endogenous'].strip() == 't']

    for d in data:

        if is_positive(d['type']) or is_positive(d['action']):

            effect = 1

        elif is_negative(d['type']) or is_negative(d['action']):

            effect = -1

        else:

            effect = 0

        for ligand_taxon in d['ligand_species'].split('|'):

            ligand_taxid = get_taxid(ligand_taxon)

            ligands = d['ligand_gene_symbol'] or d['ligand_pubchem_sid']
            ligands = ligands.split('|')
            targets = (
                d['target_uniprot'] or
                d['target_ligand_uniprot'] or
                d['target_ligand_pubchem_sid']
            )
            targets = targets.split('|')
            references = d['pubmed_id'].split('|') if d['pubmed_id'] else []

            if process_interactions:

                for ligand, target in itertools.product(ligands, targets):

                    interactions.append(
                        GuideToPharmacologyInteraction(
                            ligand = ligand,
                            ligand_id_type = (
                                'genesymbol'
                                    if d['ligand_gene_symbol'] else
                                'pubchem_sid'
                                    if d['ligand_pubchem_sid'] else
                                None
                            ),
                            target = target,
                            target_id_type = (
                                'uniprot'
                                    if (
                                        d['target_uniprot'] or
                                        d['target_ligand_uniprot']
                                    ) else
                                'pubchem_sid'
                                    if d['target_ligand_pubchem_sid'] else
                                None
                            ),
                            target_is_ligand = bool(d['target_ligand']),
                            ligand_organism = ligand_taxid,
                            target_organism = get_taxid(d['target_species']),
                            effect = effect,
                            ligand_location = (
                                d['ligand_context'].strip().lower() or None
                            ),
                            target_type = (
                                d['receptor_site'].strip().lower() or None
                            ),
                            ligand_endogenous = (
                                d['endogenous'].strip() == 't'
                            ),
                            pubmed_ids = references,
                        )
                    )

            if process_complexes:

                if (
                    len(targets) > 1 and (
                        d['target_uniprot'] or
                        d['target_ligand_uniprot']
                    )
                ):

                    cplex = intera.Complex(
                        components = targets,
                        sources = 'Guide2Pharma',
                        references = references,
                    )
                    key = cplex.__str__()

                    if key in complexes:

                        complexes[key] += cplex

                    else:

                        complexes[key] = cplex

                if (
                    len(ligands) > 1 and
                    d['ligand_gene_symbol']
                ):

                    ligand_uniprots = [
                        mapping.map_name0(ligand, 'genesymbol', 'uniprot')
                        for ligand in ligands
                    ]
                    ligand_uniprots = [u for u in ligand_uniprots if u]

                    if len(ligand_uniprots) > 1:

                        cplex = intera.Complex(
                            components = ligand_uniprots,
                            sources = 'Guide2Pharma',
                            references = references,
                        )
                        key = cplex.__str__()

                        if key in complexes:

                            complexes[key] += cplex

                        else:

                            complexes[key] = cplex

    return interactions, complexes


def guide2pharma_interactions(**kwargs):

    interactions, complexes = get_guide2pharma(
        process_complexes = False,
        **kwargs
    )

    return interactions


def guide2pharma_complexes(**kwargs):

    interactions, complexes = get_guide2pharma(
        process_interactions = False,
        **kwargs
    )

    return complexes


def cellphonedb_ligands_receptors():
    """
    Retrieves the set of ligands and receptors from CellPhoneDB.
    Returns tuple of sets.
    """

    receptors = set()
    ligands   = set()
    
    proteins = cellphonedb_protein_annotations()
    complexes = cellphonedb_complex_annotations()
    
    for _id, annot in itertools.chain(
        iteritems(proteins),
        iteritems(complexes)
    ):
        
        if annot.receptor:
            
            receptors.add(_id)
        
        if annot.secreted or (
            not annot.receptor and (
                annot.transmembrane or
                annot.peripheral
            )
        ):
            
            ligands.add(_id)
    
    return ligands, receptors


def _cellphonedb_annotations(url, name_method):


    def get_bool(rec, attr):

        return attr in rec and rec[attr] == 'True'


    def get_desc(rec, attr):

        desc = '%s_desc' % attr
        
        return (
            '' if (
                attr in rec and rec[attr] == 'False' or
                attr not in rec and not rec[desc]
            ) else
            rec[desc] if rec[desc] else
            attr.capitalize()
        )

    record = CellPhoneDBAnnotation

    annot = {}

    c = curl.Curl(url, large = True, silent = False)
    tab = list(csv.DictReader(c.result))

    for rec in tab:
        
        names = name_method(rec)
        
        if isinstance(names, (common.basestring, intera.Complex)):
            
            names = (names,)
        
        for name in names:
            
            annot[name] = record(
                receptor = get_bool(rec, 'receptor'),
                receptor_class = (
                    get_desc(rec, 'receptor').replace('_add', '').lower() or
                    None
                ),
                peripheral = get_bool(rec, 'peripheral'),
                secreted = get_bool(rec, 'secreted'),
                secreted_class = (
                    get_desc(rec, 'secreted').replace(' | ', ',').lower() or
                    None
                ),
                transmembrane = get_bool(rec, 'transmembrane'),
                integrin = get_bool(rec, 'integrin'),
            )

    return annot


def cellphonedb_protein_annotations(add_complex_annotations = True):
    """
    :arg bool add_complex_annotations:
        Copy the annotations of complexes to each of their member proteins.
    """
    
    def name_method(rec):
        
        uniprot = rec['uniprot']
        uniprot = _cellphonedb_hla(uniprot)
        uniprot = mapping.map_names(uniprot, 'uniprot', 'uniprot')
        
        return uniprot
    
    
    protein_annotations = _cellphonedb_annotations(
        url = urls.urls['cellphonedb_git']['proteins'],
        name_method = name_method,
    )
    
    if add_complex_annotations:
        
        complex_annotations = cellphonedb_complex_annotations()
        
        for cplex, a_cplex in iteritems(complex_annotations):
            
            for uniprot in cplex.components.keys():
                
                if uniprot in protein_annotations:
                    
                    protein_annotations[uniprot] = CellPhoneDBAnnotation(
                        *(
                            p or c if isinstance(p, bool) else p
                            for p, c in
                            zip(
                                protein_annotations[uniprot],
                                a_cplex
                            )
                        )
                    )
                    
                else:
                    
                    protein_annotations[uniprot].add(a_cplex)
    
    return protein_annotations


def _cellphonedb_hla(uniprot):
    """
    Returns *set*.
    """
    
    uniprots = None
    
    # for HLA genes in the uniprot column we have "HLA..." gene symbols
    # but not always in the standard format (with dash)
    if uniprot.startswith('HLA') and '-' not in uniprot:
        
        genesymbol = 'HLA-%s' % uniprot[3:]
        uniprots = mapping.map_name(genesymbol, 'genesymbol', 'uniprot')
    
    return uniprots or {uniprot}


def cellphonedb_complex_annotations():
    
    
    def get_uniprots(rec):

        return tuple(
            uniprot
            for uniprot in
            (rec['uniprot_%u' % i] for i in xrange(1, 5))
            if uniprot
        )


    def get_stoichiometry(rec):

        if not rec['stoichiometry']:

            return get_uniprots(rec)

        return tuple(
            mapping.map_name0(genesymbol, 'genesymbol', 'uniprot')
            for genesymbol in rec['stoichiometry'].split(';')
        )
    
    
    def name_method(rec):
        
        comp = get_stoichiometry(rec)

        cplex = intera.Complex(
            name = rec['complex_name'],
            components = comp,
            sources = 'CellPhoneDB',
            ids = rec['complex_name'],
        )
        
        return cplex


    return _cellphonedb_annotations(
        url = urls.urls['cellphonedb_git']['complexes'],
        name_method = name_method,
    )


def _cellphonedb_get_entity(name, complexes):
    
    if name in complexes:
        
        return (complexes[name],)
    
    if ':' in name:

        name = name.split(':')[1]

    if '_' in name:

        name = mapping.map_name0(name, 'name-entry', 'name')
        
    if not uniprot_input.is_uniprot(name):
        
        uniprot = mapping.map_name0(name, 'genesymbol', 'uniprot')
        
        name = uniprot or name

    name = _cellphonedb_hla(name)

    return (name,) if isinstance(name, common.basestring) else name


def cellphonedb_interactions():
    
    
    def get_type(entity):
        
        return (
            'ligand'
                if entity in ligands else
            'receptor'
                if entity in receptors else
            'unknown'
        )
    
    
    CellphonedbInteraction = collections.namedtuple(
        'CellphonedbInteraction',
        [
            'id_a',
            'id_b',
            'sources',
            'references',
            'interaction_type',
            'type_a',
            'type_b',
        ]
    )


    repmid = re.compile(r'PMID: ([0-9]+)')
    recomma = re.compile(r'[,;]')


    ligands, receptors = cellphonedb_ligands_receptors()
    complexes = dict(
        (
            _id,
            cplex
        )
        for cplex in cellphonedb_complexes().values()
        for _id in cplex.ids['CellPhoneDB']
    )

    url = urls.urls['cellphonedb_git']['interactions']

    c = curl.Curl(url, silent = False, large = True)

    reader = csv.DictReader(c.result)

    for rec in reader:
        
        _partner_a = _cellphonedb_get_entity(
            rec['partner_a'],
            complexes = complexes,
        )
        _partner_b = _cellphonedb_get_entity(
            rec['partner_b'],
            complexes = complexes,
        )
        
        for partner_a, partner_b in itertools.product(_partner_a, _partner_b):
            
            type_a = get_type(partner_a)
            type_b = get_type(partner_b)
            rev = partner_b == 'ligand' and partner_b == 'receptor'
            _type_a = type_b if rev else type_a
            _type_b = type_a if rev else type_b

            sources = (
                'CellPhoneDB'
                    if rec['annotation_strategy'] == 'curated' else
                '%s;CellPhoneDB' % (
                    ';'.join(
                        '%s_CellPhoneDB' % (
                            res.replace(
                                'guidetopharmacology.org',
                                'Guide2Pharma'
                            )
                        )
                        for res in
                        recomma.split(rec['annotation_strategy'])
                    )
                )
            )
            refs   = ';'.join(repmid.findall(rec['source']))
            
            type_b if rev else type_b
            
            yield(
                CellphonedbInteraction(
                    id_a = partner_b if rev else partner_a,
                    id_b = partner_a if rev else partner_b,
                    sources = sources,
                    references = refs,
                    interaction_type = '%s-%s' % (_type_a, _type_b),
                    type_a = _type_a,
                    type_b = _type_b,
                )
            )


def cellphonedb_complexes():

    annot = cellphonedb_complex_annotations()

    complexes = {}
    
    for cplex in annot.keys():
        
        key = cplex.__str__()

        if key in annot:

            cplex.add_attr('CellPhoneDB', annot[key])

        complexes[key] = cplex

    return complexes


def open_pubmed(pmid):
    """
    Opens PubMed record in web browser.

    @pmid : str or int
        PubMed ID
    """
    
    pmid = str(pmid)
    url = urls.urls['pubmed']['url'] % pmid
    webbrowser.open(url)


def only_pmids(idList, strict = True):
    """
    Return elements unchanged which comply with the PubMed ID format,
    and attempts to translate the DOIs and PMC IDs using NCBI
    E-utils.
    Returns list containing only PMIDs.

    @idList : list, str
        List of IDs or one single ID.
    @strict : bool
        Whether keep in the list those IDs which are not PMIDs,
        neither DOIs or PMC IDs or NIH manuscript IDs.
    """
    if type(idList) in common.simple_types:
        idList = [idList]
    pmids = {i for i in idList if isinstance(i, int) or i.isdigit()}
    pmcids = [i for i in idList if i.startswith('PMC')]
    dois = [i for i in idList if '/' in i]
    manuscids = [i for i in idList if i.startswith('NIHMS')]
    if not strict:
        pmids = set(pmids) | set(dois) | set(pmcids) | set(manuscids)
    if len(pmcids) > 0:
        pmids = pmids | set(pmids_list(pmcids))
    if len(dois) > 0:
        pmids = pmids | set(pmids_list(dois))
    return list(pmids)


def get_pmid(idList):
    """
    For a list of doi or PMC IDs
    fetches the corresponding PMIDs.
    """
    if type(idList) in common.simple_types:
        idList = [idList]
    url = urls.urls['pubmed-eutils']['conv'] % ','.join(str(i) for i in idList)
    c = curl.Curl(url, silent = True)
    data = c.result
    try:
        js = json.loads(data)
    except:
        js = {}
    return js


def pmids_dict(idList):
    jsn = get_pmid(idList)
    result = {'doi': {}, 'pmc': {}}
    if 'records' in jsn:
        for r in jsn['records']:
            if 'pmid' in r:
                if 'doi' in r:
                    result['doi'][r['pmid']] = r['doi']
                if 'pmcid' in r:
                    result['pmc'][r['pmid']] = r['pmcid']
    return result


def pmids_list(idList):
    jsn = get_pmid(idList)
    result = []
    if 'records' in jsn:
        for r in jsn['records']:
            if 'pmid' in r:
                result.append(r['pmid'])
    return result


def get_hprd(in_vivo = True):
    """
    Downloads and preprocesses HPRD data.
    """
    url = urls.urls['hprd_all']['url']
    files = [urls.urls['hprd_all']['ptm_file']]
    c = curl.Curl(url, silent = False, files_needed = files)
    data = c.result
    if len(data) == 0:
        return []
    data = [l.split('\t') for l in data[files[0]].split('\n')][:-1]
    if in_vivo:
        data = [i for i in data if 'in vivo' in i[9].split(';')]
    return data


def hprd_interactions(in_vivo = True):
    """
    Processes HPRD data and extracts interactions.
    Returns list of interactions.
    """
    return [i for i in get_hprd(in_vivo = in_vivo) if i[6] != '-']


def hprd_htp():

    url = urls.urls['hprd_all']['url']
    fname = urls.urls['hprd_all']['int_file']
    c = curl.Curl(url, silent = False, large = True, files_needed = [fname])

    return list(
        map(
            lambda l: l.split('\t'),
            c.result[fname].read().decode('ascii').split('\n')
        )
    )


def get_hprd_ptms(in_vivo = True):
    """
    Processes HPRD data and extracts PTMs.
    Returns list of kinase-substrate interactions.
    """
    ptms = []
    non_digit = re.compile(r'[^\d]+')
    data = get_hprd(in_vivo = in_vivo)
    for ptm in data:
        if ptm[6] != '-':
            resnums = [
                int(nn)
                for nn in [non_digit.sub('', n) for n in ptm[4].split(';')]
                if len(nn) > 0
            ]
            for resnum in resnums:
                ptms.append({
                    'resaa': ptm[5],
                    'resnum': resnum,
                    'typ': ptm[8].lower(),
                    'references': ptm[10].split(','),
                    'kinase': ptm[6],
                    'substrate_refseqp': ptm[3],
                    'substrate': ptm[1],
                    'start': max(resnum - 7, 1),
                    'end': resnum + 7,
                    'instance': None
                })
    return ptms


def disgenet_annotations(dataset = 'curated'):
    """
    Downloads and processes the list of all human disease related proteins
    from DisGeNet.
    Returns dict of dicts.

    @dataset : str
        Name of DisGeNet dataset to be obtained:
        `curated`, `literature`, `befree` or `all`.
    """

    DisGeNetAnnotation = collections.namedtuple(
        'DisGeNetAnnotation',
        [
            'disease',
            'score',
            'dsi',
            'dpi',
            'nof_pmids',
            'nof_snps',
            'source',
        ]
    )

    url = urls.urls['disgenet']['url'] % dataset
    c = curl.Curl(
        url,
        silent = False,
        large = True,
        encoding = 'utf-8',
        default_mode = 'r',
    )
    reader = csv.DictReader(c.result, delimiter = '\t')
    data = collections.defaultdict(set)

    for rec in reader:

        uniprots = mapping.map_name(
            rec['geneSymbol'],
            'genesymbol',
            'uniprot',
        )

        if not uniprots:

            continue
        
        for uniprot in uniprots:
            
            data[uniprot].add(
                DisGeNetAnnotation(
                    disease = rec['diseaseName'],
                    score = float(rec['score']),
                    dsi = float(rec['DSI']) if rec['DSI'] else None,
                    dpi = float(rec['DPI']) if rec['DPI'] else None,
                    nof_pmids = int(rec['NofPmids']),
                    nof_snps = int(rec['NofSnps']),
                    source = tuple(
                        x.strip()
                        for x in rec['source'].split(';')
                    ),
                )
            )

    return data


def load_lmpid(organism = 9606):
    """
    Reads and processes LMPID data from local file
    `pypath.data/LMPID_DATA_pubmed_ref.xml`.
    The file was provided by LMPID authors and is now
    redistributed with the module.
    Returns list of domain-motif interactions.
    """
    result = []
    
    url = urls.urls['lmpid']['url']
    c = curl.Curl(url, silent = False, large = False)
    
    soup = bs4.BeautifulSoup(c.result, 'html.parser')
    uniprots = uniprot_input.get_db(organism = organism, swissprot = None)
    prg = progress.Progress(
        len(soup.find_all('record')), 'Processing data from LMPID', 21)
    
    for rec in soup.find_all('record'):
        
        prg.step()
        uniprot_bait = rec.bait_uniprot_id.text
        uniprot_prey = rec.prey_uniprot_id.text
        
        if uniprot_bait in uniprots and uniprot_prey in uniprots:
            
            result.append({
                'bait': uniprot_bait,
                'prey': uniprot_prey,
                'refs': [x.strip() for x in rec.references.text.split(',')],
                'pos':
                [int(x) for x in rec.sequence_position.text.split('-')],
                'inst': rec.motif_instance.text,
                'dom': rec.interacting_domain.text
            })
    
    prg.terminate()
    
    return result


def lmpid_interactions(organism = 9606):
    """
    Converts list of domain-motif interactions supplied by
    `pypath.dataio.load_lmpid()` to list of interactions.
    """
    data = load_lmpid(organism = organism)
    return [[l['prey'], l['bait'], ';'.join(l['refs'])] for l in data]


def lmpid_dmi(organism = 9606):
    """
    Converts list of domain-motif interactions supplied by
    `pypath.dataio.load_lmpid()` to list of
    `pypath.intera.DomainMotif() objects.
    """
    data = load_lmpid(organism = organism)
    return [{
        'motif_protein': l['bait'],
        'domain_protein': l['prey'],
        'instance': l['inst'],
        'motif_start': l['pos'][0],
        'motif_end': l['pos'][1],
        'domain_name': l['dom'],
        'domain_name_type': 'name',
        'refs': l['refs']
    } for l in data]


def get_hsn():
    """
    Downloads and processes HumanSignalingNetwork version 6
    (published 2014 Jan by Edwin Wang).
    Returns list of interactions.
    """
    url = urls.urls['hsn']['url']
    c = curl.Curl(url, silent = False, large = True)
    data = c.result
    data = [r.split(',') for r in data if len(r) > 0]
    return data


def get_li2012():
    """
    Reads supplementary data of Li 2012 from local file.
    Returns table (list of lists).
    """
    url = urls.urls['li2012']['url']
    c = curl.Curl(url, silent = False, large = True)
    xls = c.fileobj
    xlsfile = xls.name
    xls.close()
    tbl = inputs_common.read_xls(xlsfile, sheet = 'File S1')
    return filter(lambda l: len(l[-1]) > 0, map(lambda l: l[:7], tbl[2:]))


def li2012_interactions():
    """
    Converts table read by `pypath.dataio.get_li2012()` to
    list of interactions.
    """
    result = []
    data = get_li2012()
    for l in data:
        subs_protein = l[1].split('/')[0]
        tk_protein = l[2].split()[0]
        reader_protein = l[3].split()[0]
        route = l[4]
        result.append((tk_protein, subs_protein, route, 'phosphorylation'))
        result.append((subs_protein, reader_protein, route,
                       'phosphomotif_binding'))
    return [list(l) for l in common.uniq_list(result)]


def li2012_phospho():
    """
    Converts table read by `pypath.dataio.get_li2012()` to
    list of dicts of kinase-substrate interactions.
    """
    
    result = []
    non_digit = re.compile(r'[^\d]+')
    data = get_li2012()
    
    for l in data:
        
        subs_protein = l[1].split('/')[0]
        tk_protein = l[2].split()[0]
        subs_resnum = int(non_digit.sub('', l[1].split('/')[1]))
        result.append(
            (
                subs_protein, # substrate
                tk_protein, # kinase
                None, # instance
                None, # start
                None, # end
                'Y',  # residue letter
                subs_resnum, # residue offset
            )
        )
    
    result = [
        dict(
            zip(
                [
                    'substrate', 'kinase', 'instance',
                    'start', 'end', 'resaa',
                    'resnum',
                ],
                list(l)
            )
        )
        for l in common.uniq_list(result)
    ]
    
    return result


def li2012_dmi():
    """
    Converts table read by `pypath.dataio.get_li2012()` to
    list of `pypath.intera.DomainMotif()` objects.
    Translates GeneSymbols to UniProt IDs.

    @mapper : pypath.mapping.Mapper()
        If not provided, a new `Mapper()` instance will be
        initialized, reserving more memory.
    """

    result = {}
    nondigit = re.compile(r'[^\d]+')
    se = uniprot_input.swissprot_seq(isoforms = True)
    data = get_li2012()

    for l in data:
        subs_protein = l[1].split('/')[0]
        tk_protein = l[2].split()[0]
        reader_protein = l[3].split()[0]
        subs_uniprots = mapping.map_name(subs_protein, 'genesymbol', 'uniprot')
        tk_uniprots = mapping.map_name(tk_protein, 'genesymbol', 'uniprot')
        reader_uniprots = mapping.map_name(reader_protein, 'genesymbol',
                                          'uniprot')
        subs_resnum = int(non_digit.sub('', l[1].split('/')[1]))
        for su in subs_uniprots:
            if su in se:
                subs_iso = None
                for iso, s in iteritems(se[su].isof):
                    if se[su].get(subs_resnum, isoform = iso) == 'Y':
                        subs_iso = iso
                        break
                if subs_iso:
                    start = min(1, subs_resnum - 7)
                    end = max(subs_resnum + 7, len(se[su].isof[subs_iso]))
                    for ku in tk_uniprots:
                        res = intera.Residue(
                            subs_resnum, 'Y', su, isoform = subs_iso)
                        mot = intera.Motif(
                            su,
                            start,
                            end,
                            isoform = subs_iso,
                            instance = se[su].get(start, end, isoform = subs_iso))
                        ptm = intera.Ptm(su,
                                         motif = mot,
                                         residue = res,
                                         isoform = subs_iso,
                                         source = 'Li2012')
                        dom = intera.Domain(ku)
                        dommot = intera.DomainMotif(
                            domain = dom, ptm = ptm, sources = ['Li2012'])
                        result = {}
    return result


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

    cachefile = cachefile or settings.get('trip_preprocessed')

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
        'Functional consequence': 'fc'
    }

    interactors = {}
    base_url = urls.urls['trip']['base']
    show_url = urls.urls['trip']['show']
    # url = urls.urls['trip']['url']
    # json_url = urls.urls['trip']['json']
    # c = curl.Curl(json_url, silent = False)
    # jsn = c.result
    # c = curl.Curl(url, silent = False)
    # data = c.result
    # jsn = json.loads(jsn, encoding = 'utf-8')
    c = curl.Curl(base_url)
    mainhtml = c.result
    mainsoup = bs4.BeautifulSoup(mainhtml, 'html.parser')
    trppages = common.flat_list(
        [[a.attrs['href'] for a in ul.find_all('a')]
         for ul in mainsoup.find(
             'div', id = 'trp_selector').find('ul').find_all('ul')])

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

    url = urls.urls['trip']['show'] % syn
    c = curl.Curl(url)
    html = c.result
    soup = bs4.BeautifulSoup(html, 'html.parser')

    return trip_find_uniprot(soup)


def trip_find_uniprot(soup):
    """
    Looks up a UniProt name in table downloaded from TRIP
    webpage.

    @soup : bs4.BeautifulSoup
        The `BeautifulSoup` instance returned by `pypath.dataio.trip_get_uniprot()`.
    """
    for tr in soup.find_all('div', id = 'tab2')[0].find_all('tr'):
        if tr.find('td') is not None and tr.find('td').text.strip() == 'Human':
            uniprot = tr.find_all('td')[2].text.strip()
            return uniprot
    return None


def trip_process(exclude_methods = ['Inference', 'Speculation'],
                 predictions = False,
                 species = 'Human',
                 strict = False):
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
    for uniprots in common.uniq_list(
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


def trip_interactions(exclude_methods = ['Inference', 'Speculation'],
                      predictions = False,
                      species = 'Human',
                      strict = False):
    """
    Obtains processed TRIP interactions by calling `pypath.dataio.trip_process()`
    and returns list of interactions. All arguments are passed to
    `trip_process()`, see their definition there.
    """
    data = trip_process(exclude_methods, predictions, species, strict)

    def trip_effect(eff):
        pos = set([
            'Sensitization', 'Activation', 'Increase in plasma membrane level',
            'Increase in lysosomal membrane level', 'New channel creation'
        ])
        neg = set([
            'Desensitization', 'Decrease in plasma membrane level',
            'Inhibition', 'Internalization from membrane by ligand',
            'Retain in the endoplasmic reticulum'
        ])
        return 'stimulation' if len(eff & pos) > 0 \
            else 'inhibition' if len(eff & neg) > 0 else 'unknown'

    return [[
        unipr[0], unipr[1], ';'.join(d['refs']), ';'.join(d['methods']),
        trip_effect(d['effect'])
    ] for unipr, d in iteritems(data)]


def load_macrophage():
    """
    Loads Macrophage from local file.
    Returns list of interactions.
    """
    fname = urls.files['macrophage']
    fname = os.path.join(common.ROOT, 'data', fname)
    with open(fname, 'r') as f:
        data = f.read()
    data = data.replace('?', '').replace('->', ',')


def kegg_interactions():
    """
    Downloads and processes KEGG Pathways.
    Returns list of interactions.
    """

    rehsa = re.compile(r'.*(hsa[0-9]+).*')
    req_hdrs = [
        'Referer: http://www.genome.jp/kegg-bin/show_pathway'
        '?map=hsa04710&show_description=show'
    ]
    hsa_list = []
    interactions = []

    c = curl.Curl(urls.urls['kegg_pws']['list_url'], silent = True)
    htmllst = c.result
    lstsoup = bs4.BeautifulSoup(htmllst, 'html.parser')

    for a in lstsoup.find_all('a', href = True):
        m = rehsa.match(a['href'])
        if m:
            hsa_list.append((m.groups(0)[0], a.text))

    prg = progress.Progress(
        len(hsa_list), 'Processing KEGG Pathways', 1, percent = False)

    for hsa, pw in hsa_list:

        prg.step()
        c = curl.Curl(urls.urls['kegg_pws']['kgml_url_2'] % hsa,
                      silent = True,
                      req_headers = req_hdrs)
        kgml = c.result
        kgmlsoup = bs4.BeautifulSoup(kgml, 'html.parser')
        entries = {}

        for ent in kgmlsoup.find_all('entry'):
            gr = ent.find('graphics')
            if gr and 'name' in gr.attrs:
                entries[ent.attrs['id']] = [
                    n.strip()
                    for n in gr.attrs['name'].replace('...', '').split(',')
                ]

        uentries = dict([(eid, common.uniq_list(
            common.flat_list([
                mapping.map_name(
                    gn, 'genesymbol', 'uniprot', strict = True) for gn in gns
            ]))) for eid, gns in iteritems(entries)])

        for rel in kgmlsoup.find_all('relation'):
            st = rel.find('subtype')
            if (
                rel.attrs['entry1'] in uentries and
                rel.attrs['entry2'] in uentries and
                st and
                'name' in st.attrs
            ):

                for u1 in uentries[rel.attrs['entry1']]:

                    for u2 in uentries[rel.attrs['entry2']]:

                        interactions.append((u1, u2, st.attrs['name'], pw))

    prg.terminate()
    return common.uniq_list(interactions)


def kegg_pathways():

    data = kegg_interactions()
    pws = common.uniq_list(map(lambda i: i[3], data))
    proteins_pws = dict(map(lambda pw: (pw, set([])), pws))
    interactions_pws = dict(map(lambda pw: (pw, set([])), pws))
    
    for u1, u2, eff, pw in data:
        
        proteins_pws[pw].add(u1)
        proteins_pws[pw].add(u2)
        interactions_pws[pw].add((u1, u2))
    
    return proteins_pws, interactions_pws


def kegg_pathway_annotations():

    KeggPathway = collections.namedtuple(
        'KeggPathway', ['pathway'],
    )


    result = collections.defaultdict(set)

    proteins, interactions = kegg_pathways()

    for pathway, uniprots in iteritems(proteins):

        record = KeggPathway(pathway = pathway)

        for uniprot in uniprots:

            result[uniprot].add(record)

    return result


def kinasedotcom_annotations():
    """
    Downloads and processes kinase annotations from kinase.com.
    """
    
    KinasedotcomAnnotation = collections.namedtuple(
        'KinasedotcomAnnotation',
        ['group', 'family', 'subfamily']
    )
    KinasedotcomAnnotation.__new__.__defaults__ = (None,)
    
    
    def add_record(uniprot, rec, offset = 2):
        
        if rec[offset].strip():
            
            result[uniprot].add(
                KinasedotcomAnnotation(
                    group = rec[offset].strip(),
                    family = rec[offset + 1].strip(),
                    subfamily = rec[offset + 2].strip() or None,
                )
            )
    
    
    url = urls.urls['kinome']['url']
    c = curl.Curl(url, large = True, silent = False)
    xlsf = c.fileobj
    xlsname = xlsf.name
    xlsf.close()
    tbl = inputs_common.read_xls(xlsname)
    
    result = collections.defaultdict(set)
    
    for rec in tbl:
        
        uniprots = mapping.map_name(rec[23].strip(), 'genesymbol', 'uniprot')
        
        for uniprot in uniprots:
            
            add_record(uniprot, rec)
            
            if rec[12].strip():
                
                add_record(uniprot, rec, offset = 12)
    
    return result


def phosphatome_annotations():
    """
    Downloads the list of phosphatases from Chen et al, Science Signaling
    (2017) Table S1.
    """

    PhosphatomeAnnotation = collections.namedtuple(
        'PhosphatomeAnnotation',
        [
            'fold',
            'family',
            'subfamily',
            'has_protein_substrates',
            'has_non_protein_substrates',
            'has_catalytic_activity',
        ],
    )

    url = urls.urls['phosphatome']['url']
    c = curl.Curl(url, large = True, silent = False, default_mode = 'rb')
    tbl = inputs_common.read_xls(c.result['aag1796_Tables S1 to S23.xlsx'])

    data = collections.defaultdict(set)

    for rec in tbl[2:]:

        uniprots = mapping.map_name(rec[0], 'genesymbol', 'uniprot')
        
        for uniprot in uniprots:
            
            data[uniprot].add(
                PhosphatomeAnnotation(
                    fold = rec[2],
                    family = rec[3],
                    subfamily = rec[4],
                    has_protein_substrates = rec[21].strip().lower() == 'yes',
                    has_non_protein_substrates = (
                        rec[22].strip().lower() == 'yes'
                    ),
                    has_catalytic_activity = rec[23].strip().lower() == 'yes',
                )
            )

    return data


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


def dgidb_annotations():
    """
    Downloads druggable protein annotations from DGIdb.
    """
    
    DgidbAnnotation = collections.namedtuple(
        'DgidbAnnotation',
        ['category'],
    )
    
    
    url = urls.urls['dgidb']['categories']
    
    c = curl.Curl(url = url, silent = False, large = True)
    
    data = csv.DictReader(c.result, delimiter = '\t')
    
    result = collections.defaultdict(set)
    
    for rec in data:
        
        uniprots = mapping.map_name(
            rec['entrez_gene_symbol'],
            'genesymbol',
            'uniprot',
        )
        
        for uniprot in uniprots:
            
            result[uniprot].add(
                DgidbAnnotation(
                    category = rec['category']
                )
            )
    
    return result


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
        url = urls.urls['reactome']['biopax_l3']
        c = curl.Curl(
            url,
            silent = False,
            large = True,
            files_needed = ['%s.owl' % organisms[organism]]).values()[0]
        with open(unzipped, 'w') as _unzipped:
            while True:
                chunk = c.result.read(4096)
                if not chunk:
                    break
                _unzipped.write(chunk)
        c.result.close()
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
        return pickle.load(open(cachefile, 'r'))
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
    if type(biopax_file) is file:
        bp_filesize = os.path.getsize(biopax_file.name)
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
        dict(controls.items() + catalyses.items()), entity_uniprot,
        dict(reactions_uniprots.items() + complexassemblies_uniprots.items()),
        publications)
    for caref, ca in iteritems(complexassemblies_uniprots):
        controls_uniprots[caref] = {
            'type': 'BINDING',
            'refs':
            [publications[r] for r in ca['refs'] if r in publications],
            'controller': None,
            'controlled': ca
        }
    del entity_uniprot
    pickle.dump(controls_uniprots, open(cachefile, 'w'))
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


def reactome_interactions(cacheFile = None, **kwargs):
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
    else:
        while True:
            sys.stdout.write(
                '\nProcessing Reactome requires huge memory.\n'
                'Please hit `y` if you have at least 2G free memory,\n'
                'or `n` to omit Reactome.\n'
                'After processing once, it will be saved in \n'
                '%s, so next time can be loaded quickly.\n\n'
                'Process Reactome now? [y/n]\n' % cacheFile)
            sys.stdout.flush()
            answer = raw_input().lower()
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

def signalink_interactions():
    """
    Reads and processes SignaLink3 interactions from local file.
    Returns list of interactions.
    """

    repar = re.compile(r'.*\(([a-z\s]+)\)')
    repref = re.compile(r'(?:.*:)?((?:[\w]+[^\s])?)\s?')
    notNeeded = set(['acsn', 'reactome'])
    nodes = {}
    interactions = []

    def _get_attr(attrs, attrName):
        return _process_attr(attrs[attrName]) if attrName in attrs else ''

    def _process_attr(attr):
        m = repar.match(attr)
        if m is not None:
            return m.groups()[0]
        else:
            return attr
    
    url_nodes = urls.urls['signalink']['nodes']
    c_nodes = curl.Curl(url_nodes, silent = False, large = True)
    url_edges = urls.urls['signalink']['edges']
    c_edges = curl.Curl(url_edges, silent = False, large = True)
    
    for l in c_nodes.result:

        if len(l) > 0:

            l = l.split('\t')
            _id = int(l[0])
            uniprot = repref.sub('\\1', l[1])
            pathways = [
                pw.split(':')[-1].strip() for pw in l[4].split('|')
                if pw.split(':')[0] not in notNeeded
            ]
            nodes[_id] = [uniprot, pathways]

    lPrev = None

    for l in c_edges.result:

        l = l.strip().split('\t')

        if lPrev is not None:
            l = lPrev + l[1:]
            lPrev = None

        if len(l) == 13:

            if l[-1] == '0':

                dbs = [
                    _process_attr(db.split(':')[-1])
                    for db in l[9].replace('"', '').split('|')
                ]
                dbs = list(set(dbs) - notNeeded)
                if len(dbs) == 0:
                    continue
                idSrc = int(l[1])
                idTgt = int(l[2])

                uniprotSrc = repref.sub('\\1', l[3])
                uniprotTgt = repref.sub('\\1', l[4])

                if not uniprotSrc or not uniprotTgt:

                    continue

                refs = [ref.split(':')[-1] for ref in l[7].split('|')]
                attrs = dict(
                    tuple(attr.strip().split(':', 1))
                    for attr in l[8].replace('"', '').split('|'))
                interactions.append([
                    uniprotSrc, uniprotTgt, ';'.join(refs), ';'.join(dbs),
                    _get_attr(attrs, 'effect'),
                    _get_attr(attrs, 'is_direct'),
                    _get_attr(attrs, 'is_directed'),
                    _get_attr(attrs, 'molecular_background'),
                    ';'.join(nodes[idSrc][1]), ';'.join(nodes[idTgt][1])
                ])
        else:
            lPrev = l
    
    return interactions


def signalink_pathway_annotations():

    SignalinkPathway = collections.namedtuple(
        'SignalinkPathway',
        ['pathway', 'core'],
    )


    result = collections.defaultdict(set)

    interactions = signalink_interactions()

    for i in interactions:
        
        for idx in (0, 1):
            
            for pathway in i[idx + 8].split(';'):
                
                core = 'non-core' not in pathway
                pathway = (
                    pathway.split('(')[0].strip().replace('/Wingless', '')
                )
                
                for uniprot in mapping.map_name(i[idx], 'uniprot', 'uniprot'):
                    
                    result[uniprot].add(
                        SignalinkPathway(pathway = pathway, core = core)
                    )

    return result


def _netbiol_interactions(database):
    
    url = urls.urls[database]['url']
    c = curl.Curl(url, silent = True, large = False)
    
    return [row.split(',') for row in c.result.split('\n')]


def arn_interactions():
    
    return _netbiol_interactions(database = 'arn')


def nrf2ome_interactions():
    
    return _netbiol_interactions(database = 'nrf2ome')


def get_laudanna_directions():
    """
    Downloads and processes the SignalingFlow edge attributes
    from Laudanna Lab.
    Returns list of directions.
    """
    url = urls.urls['laudanna']['sigflow_rescued']
    c = curl.Curl(url, silent = False)
    data = c.result
    data = data.split('\n')[1:]
    directions = []
    for l in data:
        if len(l) > 0:
            directions.append(l.split('=')[0].strip().split(' (pp) '))
    return directions


def get_laudanna_effects():
    """
    Downloads and processes the SignalingDirection edge attributes
    from Laudanna Lab.
    Returns list of effects.
    """
    url = urls.urls['laudanna']['sigdir_rescued']
    c = curl.Curl(url, silent = False)
    data = c.result
    data = data.split('\n')[1:]
    effects = []
    for l in data:
        if len(l) > 0:
            l = l.split('=')
            effects.append(l[0].strip().split(' (pp) ') + [l[1].strip()])
    return effects


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


def wang_interactions():
    """
    Downloads and processes Wang Lab HumanSignalingNetwork.
    Returns list of interactions.
    """

    url = urls.urls['wang']['url']
    c = curl.Curl(url, silent = False)
    data = c.result
    data = data.split('\n')
    effects = []
    nodes = {}
    reading_nodes = False
    reading_edges = False

    for l in data:

        if len(l.strip()) == 0:
            reading_nodes = False
            reading_edges = False
        l = l.split(',')
        if reading_nodes:
            nodes[l[0]] = l[1]
        if reading_edges:
            effects.append([nodes[l[0]], nodes[l[1]], l[2]])
        if l[0].startswith('Node'):
            reading_nodes = True
        if l[0].startswith('From'):
            reading_nodes = False
            reading_edges = True

    return effects


def biogrid_interactions(organism = 9606, htp_limit = 1, ltp = True):
    """
    Downloads and processes BioGRID interactions.
    Keeps only the "low throughput" interactions.
    Returns list of interactions.

    @organism : int
        NCBI Taxonomy ID of organism.
    @htp_limit : int
        Exclude interactions only from references
        cited at more than this number of interactions.
    """
    organism = str(organism)
    interactions = []
    refc = []
    url = urls.urls['biogrid']['url']
    c = curl.Curl(url, silent = False, large = True)
    f = next(iter(c.result.values()))
    nul = f.readline()

    for l in f:

        l = l.split('\t')
        if len(l) > 17:
            if l[17].startswith('Low') or not ltp and l[15] == organism and l[
                    16] == organism:
                interactions.append([l[7], l[8], l[14]])
                refc.append(l[14])

    refc = Counter(refc)
    if htp_limit is not None:
        interactions = [i for i in interactions if refc[i[2]] <= htp_limit]
    return interactions


def acsn_interactions(keep_in_complex_interactions = True):
    """
    Processes ACSN data from local file.
    Returns list of interactions.

    @keep_in_complex_interactions : bool
        Whether to include interactions from complex expansion.
    """
    
    names_url = urls.urls['acsn']['names']
    ppi_url = urls.urls['acsn']['ppi']
    names_c = curl.Curl(names_url, silent = False, large = True)
    ppi_c = curl.Curl(ppi_url, silent = False, large = True)
    
    names = {}
    interactions = []
    
    for l in names_c.result:
        l = l.strip().split('\t')
        names[l[0]] = l[2:]
    
    _ = next(ppi_c.result)
    for l in ppi_c.result:
        l = l.strip().split('\t')
        if l[0] in names:
            for a in names[l[0]]:
                if l[2] in names:
                    for b in names[l[2]]:
                        if keep_in_complex_interactions:
                            if 'PROTEIN_INTERACTION' in l[1]:
                                l[1].replace('COMPLEX_EXPANSION',
                                             'IN_COMPLEX_INTERACTION')
                        interactions.append([a, b, l[1], l[3]])
    
    return interactions


def get_graphviz_attrs():
    """
    Downloads graphviz attribute list from graphviz.org.
    Returns 3 dicts of dicts: graph_attrs, vertex_attrs and edge_attrs.
    """
    url = urls.urls['graphviz']['url']
    c = curl.Curl(url)
    html = c.result
    soup = bs4.BeautifulSoup(html, 'lxml')
    vertex_attrs = {}
    edge_attrs = {}
    graph_attrs = {}
    for tbl in soup.find_all('table'):
        if tbl.find('tr').text.startswith('Name'):
            for r in tbl.find_all('tr'):
                r = r.find_all('td')
                if len(r) > 0:
                    usedby = r[1].text
                    this_attr = {
                        'type': r[2].text.strip(),
                        'default': r[3].text.strip(),
                        'min': r[4].text.strip(),
                        'notes': r[5].text.strip()
                    }
                    attr_name = r[0].text.strip()
                    if 'N' in usedby:
                        vertex_attrs[attr_name] = this_attr
                    if 'E' in usedby:
                        edge_attrs[attr_name] = this_attr
                    if 'G' in usedby:
                        graph_attrs[attr_name] = this_attr
            break
    return graph_attrs, vertex_attrs, edge_attrs


def get_phosphosite(cache = True, ncbi_tax_id = 9606):
    """
    Downloads curated and HTP data from Phosphosite,
    from preprocessed cache file if available.
    Processes BioPAX format.
    Returns list of interactions.
    """
    
    curated_cache = urls.files['phosphosite']['curated']
    noref_cache = urls.files['phosphosite']['noref']
    
    if cache and os.path.exists(curated_cache) and os.path.exists(noref_cache):
        
        return (
            pickle.load(open(curated_cache, 'rb')),
            pickle.load(open(noref_cache, 'rb')),
        )
    
    result_curated = []
    result_noref = []
    url = urls.urls['psite_bp']['url']
    c = curl.Curl(url, silent = False, large = True)
    bpax = c.gzfile
    xml = ET.parse(bpax)
    xmlroot = xml.getroot()
    bpprefix = '{http://www.biopax.org/release/biopax-level3.owl#}'
    rdfprefix = '{http://www.w3.org/1999/02/22-rdf-syntax-ns#}'
    proteins = {}
    for p in xmlroot.iter(bpprefix + 'ProteinReference'):
        psid = p.attrib[rdfprefix + 'ID']
        db = p.find(bpprefix + 'xref').find(bpprefix + 'UnificationXref').find(
            bpprefix + 'db').text
        up = p.find(bpprefix + 'xref').find(bpprefix + 'UnificationXref').find(
            bpprefix + 'id').text
        tax = ''
        if p.find(bpprefix + 'organism') is not None:
            tmp = p.find(bpprefix + 'organism')
            if rdfprefix + 'resource' in tmp.attrib:
                tax = tmp.attrib[rdfprefix + 'resource'].split('_')[1]
        if db == 'UniProtKB':
            up = up[0:6]
        proteins[psid] = {'id': up, 'db': db, 'species': tax, 'psid': psid}
    evidences = {}
    for p in xmlroot.iter(bpprefix + 'EvidenceCodeVocabulary'):
        evid = p.attrib[rdfprefix + 'ID'].split('_')[1]
        evname = p.find(bpprefix + 'term').text
        evidences[evid] = evname
    ev_short = {'0113': 'WB', '0427': 'MS', '0074': 'MA', '0421': 'AB'}
    nosrc = []
    notgt = []
    norefs = []
    noev = []
    noth = []
    edges = []
    
    for c in xmlroot.findall(bpprefix + 'Catalysis'):
        if rdfprefix + 'resource' in c.find(bpprefix + 'controller').attrib:
            src = 'po_' + \
                c.find(
                    bpprefix + 'controller').attrib[rdfprefix + 'resource'].split('_')[1]
        else:
            srcProt = c.find(bpprefix + 'controller').find(bpprefix +
                                                           'Protein')
            if srcProt is not None:
                src = 'po_' + srcProt.attrib[rdfprefix + 'ID'].split('_')[1]
            else:
                nosrc.append(c)
        tgtProt = c.find(bpprefix + 'controlled').iter(bpprefix +
                                                       'ProteinReference')
        tgt = next(tgtProt, None)
        if tgt is not None:
            tgt = tgt.attrib[rdfprefix + 'ID']
        else:
            tgtProt = c.find(bpprefix + 'controlled').iter(bpprefix +
                                                           'entityReference')
            tgt = next(tgtProt, None)
            if tgt is not None:
                if rdfprefix + 'resource' in tgt.attrib:
                    tgt = tgt.attrib[rdfprefix + 'resource'][1:]
            else:
                tgtProt = c.find(bpprefix + 'controlled').iter(bpprefix +
                                                               'left')
                tgt = next(tgtProt, None)
                if tgt is not None:
                    if rdfprefix + 'resource' in tgt.attrib:
                        tgt = 'po_' + \
                            tgt.attrib[rdfprefix + 'resource'].split('_')[1]
                else:
                    notgt.append(c)
        refs = c.iter(bpprefix + 'PublicationXref')
        pmids = []
        for r in refs:
            pm = r.attrib[rdfprefix + 'ID'].split('_')
            if pm[0] == 'pmid':
                pmids.append(pm[1])
        refs = c.iter(bpprefix + 'evidence')
        for r in refs:
            rrefs = r.iter(bpprefix + 'xref')
            for rr in rrefs:
                if rdfprefix + 'resource' in rr.attrib:
                    pm = rr.attrib[rdfprefix + 'resource'].split('_')
                    if pm[0] == 'pubmed':
                        pmids.append(pm[1])
        evs = []
        for e in c.iter(bpprefix + 'evidenceCode'):
            if rdfprefix + 'resource' in e.attrib:
                evs.append(ev_short[e.attrib[rdfprefix + 'resource'].split('_')
                                    [1]])
            else:
                ev = e.find(bpprefix + 'EvidenceCodeVocabulary')
                evs.append(ev_short[ev.attrib[rdfprefix + 'ID'].split('_')[1]])
        for e in c.iter(bpprefix + 'evidence'):
            if rdfprefix + 'resource' in e.attrib:
                ev = e.attrib[rdfprefix + 'resource'].split('_')
                if len(ev) == 4:
                    if len(ev[3]) == 4:
                        evs.append(ev_short[ev[3]])
        if (src is not None and tgt is not None and src in proteins and
                tgt in proteins and proteins[src]['id'] is not None and
                proteins[tgt]['id'] is not None):
            edges.append({
                'src': proteins[src],
                'tgt': proteins[tgt],
                'pmids': list(set(pmids)),
                'evs': list(set(evs))
            })
            if len(evs) == 0:
                noev.append(c)
            if len(pmids) == 0:
                norefs.append(c)
            if len(evs) == 0 and len(pmids) == 0:
                noth.append(c)
    
    if ncbi_tax_id:
        
        all_uniprots = uniprot_input.all_uniprots(organism = ncbi_tax_id)
    
    for e in edges:
        
        if (
            ncbi_tax_id and (
                e['src']['id'] not in all_uniprots or
                e['tgt']['id'] not in all_uniprots
            )
        ):
            
            continue
        
        this_iaction = [
            e['src']['id'], e['tgt']['id'], e['src']['species'],
            e['tgt']['species'], ';'.join(e['evs']), ';'.join(e['pmids'])
        ]
        
        if len(this_iaction[-1]) > 0:
            
            result_curated.append(this_iaction)
            
        else:
            
            result_noref.append(this_iaction)
    
    pickle.dump(result_curated, open(curated_cache, 'wb'))
    pickle.dump(result_noref, open(noref_cache, 'wb'))
    return result_curated, result_noref


def get_phosphosite_new(cache = True):
    """
    Downloads curated and HTP data from Phosphosite,
    from preprocessed cache file if available.
    Processes BioPAX format.
    Returns list of interactions.
    """
    
    curated_cache = urls.files['phosphosite']['curated']
    noref_cache = urls.files['phosphosite']['noref']
    
    if (
        cache and
        os.path.exists(curated_cache) and
        os.path.exists(noref_cache)
    ):
        
        with pen(curated_cache, 'rb') as fp:
            
            data_curated = pickle.load(fp)
        
        with pen(noref_cache, 'rb') as fp:
            
            data_noref = pickle.load(fp)
            
        return data_curated, data_noref
    
    
    def collect_items(tagname, process_method):
        
        result = {}
        
        for p in xmlroot.iter(tagname):
            
            key, value = process_method(p)
            result[key] = value
        
        return result
    
    
    def process_protein(protein):
        
        protein_id = protein.attrib['%sID' % rdfprefix]
        database = (
            protein.find(
                '%sxref' % bpprefix
            ).find(
                '%sUnificationXref' % bpprefix
            ).find(
                '%sdb' % bpprefix
            ).text
        )
        identifier = (
            protein.find(
                '%sxref' % bpprefix
            ).find(
                '%sUnificationXref' % bpprefix
            ).find(
                '%sid' % bpprefix
            ).text
        )
        
        organism = None
        e_organism = protein.find('%sorganism' % bpprefix)
        if (
            e_organism is not None and
            '%sresource' % rdfprefix in e_organism.attrib
        ):
            organism = (
                e_organism.attrib['%sresource' % rdfprefix].split('_')[1]
            )
        
        return protein_id, (databas, identifier, organism)
    
    def process_site(site):
        
        site_id = site.attrib['%sID' % rdfprefix]
        site_offset = site.find('%ssequencePosition').text
        
        return site_id, site_offset
    
    
    def process_modification(seqmodvoc):
        
        mod_id = seqmodvoc.attrib['%sID' % rdfprefix]
        residue, mod = mod_id.split('_').split('-')
        
        return mod_id, (residue, mod)
    
    
    def get_resource(elem, resource_tag):
        
        res_attr = '%sresource' % rdfprefix
        
        if res_attr in elem.attrib:
            
            return elem.attrib[res_attr][1:]
            
        else:
            
            return elem.find(resource_tag).attrib['%sID' % rdfprefix]
    
    
    def process_feature(feature):
        
        feature_id = feature.attrib['%sID' % rdfprefix]
        site = get_resource(
            feature.find('%sfeatureLocation' % bpprefix),
            '%sSequenceSite' % bpprefix,
        )
        modification = get_resource(
            feature.find('%smodificationType' % bpprefix),
            '%sSequenceModificationVocabulary' % bpprefix,
        )
        
        return feature_id, (site, modification)
    
    
    result_curated = []
    result_noref = []
    bpprefix = '{http://www.biopax.org/release/biopax-level3.owl#}'
    rdfprefix = '{http://www.w3.org/1999/02/22-rdf-syntax-ns#}'
    
    url = urls.urls['psite_bp']['url']
    c = curl.Curl(url, silent = False, large = True)
    bpax = c.gzfile
    xml = ET.parse(bpax)
    xmlroot = xml.getroot()
    
    
    proteins = collect_items(
        '%sProtein' % bpprefix,
        process_method = process_protein,
    )
    sites = collect_items(
        '%sSequenceSite' % bpprefix,
        process_method = process_site,
    )
    modifications = collect_items(
        '%sSequenceModificationVocabulary' % bpprefix,
        process_method = process_modification,
    )
    features = collect_items(
        '%sModificationFeature' % bpprefix,
        process_method = process_feature,
    )
    
    evidences = {}
    for p in xmlroot.iter(bpprefix + 'EvidenceCodeVocabulary'):
        evid = p.attrib[rdfprefix + 'ID'].split('_')[1]
        evname = p.find(bpprefix + 'term').text
        evidences[evid] = evname
    ev_short = {'0113': 'WB', '0427': 'MS', '0074': 'MA', '0421': 'AB'}
    nosrc = []
    notgt = []
    norefs = []
    noev = []
    noth = []
    edges = []
    
    for c in xmlroot.findall(bpprefix + 'Catalysis'):
        if rdfprefix + 'resource' in c.find(bpprefix + 'controller').attrib:
            src = 'po_' + \
                c.find(
                    bpprefix + 'controller').attrib[rdfprefix + 'resource'].split('_')[1]
        else:
            srcProt = c.find(bpprefix + 'controller').find(bpprefix +
                                                           'Protein')
            if srcProt is not None:
                src = 'po_' + srcProt.attrib[rdfprefix + 'ID'].split('_')[1]
            else:
                nosrc.append(c)
        tgtProt = c.find(bpprefix + 'controlled').iter(bpprefix +
                                                       'ProteinReference')
        tgt = next(tgtProt, None)
        if tgt is not None:
            tgt = tgt.attrib[rdfprefix + 'ID']
        else:
            tgtProt = c.find(bpprefix + 'controlled').iter(bpprefix +
                                                           'entityReference')
            tgt = next(tgtProt, None)
            if tgt is not None:
                if rdfprefix + 'resource' in tgt.attrib:
                    tgt = tgt.attrib[rdfprefix + 'resource'][1:]
            else:
                tgtProt = c.find(bpprefix + 'controlled').iter(bpprefix +
                                                               'left')
                tgt = next(tgtProt, None)
                if tgt is not None:
                    if rdfprefix + 'resource' in tgt.attrib:
                        tgt = 'po_' + \
                            tgt.attrib[rdfprefix + 'resource'].split('_')[1]
                else:
                    notgt.append(c)
        refs = c.iter(bpprefix + 'PublicationXref')
        pmids = []
        for r in refs:
            pm = r.attrib[rdfprefix + 'ID'].split('_')
            if pm[0] == 'pmid':
                pmids.append(pm[1])
        refs = c.iter(bpprefix + 'evidence')
        for r in refs:
            rrefs = r.iter(bpprefix + 'xref')
            for rr in rrefs:
                if rdfprefix + 'resource' in rr.attrib:
                    pm = rr.attrib[rdfprefix + 'resource'].split('_')
                    if pm[0] == 'pubmed':
                        pmids.append(pm[1])
        evs = []
        for e in c.iter(bpprefix + 'evidenceCode'):
            if rdfprefix + 'resource' in e.attrib:
                evs.append(ev_short[e.attrib[rdfprefix + 'resource'].split('_')
                                    [1]])
            else:
                ev = e.find(bpprefix + 'EvidenceCodeVocabulary')
                evs.append(ev_short[ev.attrib[rdfprefix + 'ID'].split('_')[1]])
        for e in c.iter(bpprefix + 'evidence'):
            if rdfprefix + 'resource' in e.attrib:
                ev = e.attrib[rdfprefix + 'resource'].split('_')
                if len(ev) == 4:
                    if len(ev[3]) == 4:
                        evs.append(ev_short[ev[3]])
        if (src is not None and tgt is not None and src in proteins and
                tgt in proteins and proteins[src]['id'] is not None and
                proteins[tgt]['id'] is not None):
            edges.append({
                'src': proteins[src],
                'tgt': proteins[tgt],
                'pmids': list(set(pmids)),
                'evs': list(set(evs))
            })
            if len(evs) == 0:
                noev.append(c)
            if len(pmids) == 0:
                norefs.append(c)
            if len(evs) == 0 and len(pmids) == 0:
                noth.append(c)
    for e in edges:
        this_iaction = [
            e['src']['id'], e['tgt']['id'], e['src']['species'],
            e['tgt']['species'], ';'.join(e['evs']), ';'.join(e['pmids'])
        ]
        if len(this_iaction[-1]) > 0:
            result_curated.append(this_iaction)
        else:
            result_noref.append(this_iaction)
    pickle.dump(result_curated, open(curated_cache, 'wb'))
    pickle.dump(result_noref, open(noref_cache, 'wb'))
    return result_curated, result_noref


def _phosphosite_filter_organism(psite_data, ncbi_tax_id = 9606):
    
    all_uniprots = uniprot_input.all_uniprots(organism = ncbi_tax_id)
    
    return [
        rec
        for rec in psite_data
        if rec[0] in all_uniprots and rec[1] in all_uniprots
    ]


def get_phosphosite_curated(ncbi_tax_id = 9606):
    """
    Loads literature curated PhosphoSite data,
    from preprocessed cache file if available.
    Returns list of interactions.
    """
    
    curated_cache = urls.files['phosphosite']['curated']
    if not os.path.exists(curated_cache):
        curated, noref = get_phosphosite(ncbi_tax_id = ncbi_tax_id)
        result = curated
    else:
        result = pickle.load(open(curated_cache, 'rb'))
    
    return _phosphosite_filter_organism(result, ncbi_tax_id)


def get_phosphosite_noref(ncbi_tax_id = 9606):
    """
    Loads HTP PhosphoSite data,
    from preprocessed cache file if available.
    Returns list of interactions.
    """
    
    noref_cache = urls.files['phosphosite']['noref']
    if not os.path.exists(noref_cache):
        curated, noref = get_phosphosite(ncbi_tax_id = ncbi_tax_id)
        result = noref
    else:
        result = pickle.load(open(noref_cache, 'rb'))
    
    return _phosphosite_filter_organism(result, ncbi_tax_id)


def phosphosite_directions(organism = 'human'):
    """
    From curated and HTP PhosphoSite data generates a
    list of directions.
    """
    curated, noref = get_phosphosite()
    return [
        i[:2] for i in curated + noref if i[2] == organism and i[3] == organism
    ]


def phosphosite_interactions_all():
    
    return get_phosphosite_curated() + get_phosphosite_noref()


def get_ca1():
    """
    Downloads and processes the CA1 signaling network (Ma\'ayan 2005).
    Returns list of interactions.
    """
    url = urls.urls['ca1']['url']
    c = curl.Curl(url, silent = False, files_needed = ['S1.txt'])
    data = c.result
    return filter(lambda l: len(l) == 13,
                  map(lambda l: l.strip().split(),
                      data['S1.txt'].split('\n')[1:]))


def get_ccmap(organism = 9606):
    """
    Downloads and processes CancerCellMap.
    Returns list of interactions.

    @organism : int
        NCBI Taxonomy ID to match column #7 in nodes file.
    """
    organism = '%u' % organism
    interactions = []
    nodes_url = urls.urls['ccmap']['nodes']
    edges_url = urls.urls['ccmap']['edges']
    c = curl.Curl(
        nodes_url, silent = False,
        files_needed = ['cell-map-node-attributes.txt'])
    nodes = c.result
    c = curl.Curl(
        edges_url, silent = False,
        files_needed = ['cell-map-edge-attributes.txt'])
    edges = c.result
    nodes = dict(
        map(lambda l: (l[1], l[2].split(':')),
            filter(lambda l: l[5] == 'protein' and l[6] == organism,
                   filter(lambda l: len(l) == 7,
                          map(lambda l: l.strip().split('\t'), nodes[
                              'cell-map-node-attributes.txt'].split('\n')[
                                  1:])))))
    edges = filter(lambda l: len(l) == 7,
                   map(lambda l: l.strip().split('\t'),
                       edges['cell-map-edge-attributes.txt'].split('\n')[1:]))
    for e in edges:
        if e[1] != 'IN_SAME_COMPONENT' and e[3] in nodes and e[4] in nodes:
            for src in nodes[e[3]]:
                for tgt in nodes[e[4]]:
                    interactions.append([
                        src, tgt, 'directed' if e[1] == 'STATE_CHANGE' else
                        'undirected', e[6].strip(';').replace('PUBMED:', '')
                    ])
    return interactions


def cancer_gene_census_annotations(
        user = None,
        passwd = None,
        credentials_fname = 'cosmic_credentials',
    ):
    """
    Retrieves a list of cancer driver genes (Cancer Gene Census) from
    the Sanger COSMIC (Catalogue of Somatic Mutations in Cancer) database.

    Does not work at the moment (signature does not match error).
    """

    if not user or not passwd:

        credentials = settings.get('cosmic_credentials')

        if not credentials:
            
            if os.path.exists(credentials_fname):
                
                _log(
                    'Reading COSMIC credentials '
                    'from file `%s`.' % credentials_fname
                )
                
                with open(credentials_fname, 'r') as fp:
                    
                    credentials = dict(
                        zip(
                            ('user', 'passwd'),
                            fp.read().split('\n')[:2],
                        )
                    )
            
        else:
            
            _log('COSMIC credentials provided by `settings`.')
        
        if not credentials or {'user', 'passwd'} - set(credentials.keys()):
            
            _log(
                'No credentials available for the COSMIC website. '
                'Either set the `cosmic_credentials` key in the `settings` '
                'module (e.g. `{\'user\': \'myuser\', '
                '\'passwd\': \'mypassword\'}`), or pass them directly to the '
                '`dataio.cancer_gene_census_annotations` method.'
            )
            return {}

    else:

        credentials = {'user': user, 'passwd': passwd}

    CancerGeneCensusAnnotation = collections.namedtuple(
        'CancerGeneCensusAnnotation',
        [
            'tier',
            'hallmark',
            'somatic',
            'germline',
            'tumour_types_somatic',
            'tumour_types_germline',
            'cancer_syndrome',
            'tissue_type',
            'genetics',
            'role',
            'mutation_type',
        ]
    )


    def multi_field(content):

        return (
            tuple(sorted(i.strip() for i in content.split(',')))
                if content.strip() else
            ()
        )


    url = urls.urls['cgc']['url_new']

    auth_str = base64.b64encode(
        ('%s:%s\n' % (credentials['user'], credentials['passwd'])).encode()
    )

    req_hdrs = ['Authorization: Basic %s' % auth_str.decode()]

    c = curl.Curl(
        url,
        large = False,
        silent = False,
        req_headers = req_hdrs,
        cache = False,
    )

    access_url = json.loads(c.result)

    if 'url' not in access_url:

        _log(
            'Could not retrieve COSMIC access URL. '
            'Most likely the authentication failed. '
            'The reply was: `%s`' % c.result
        )

        return None

    c = curl.Curl(
        access_url['url'],
        large = True,
        silent = False,
        bypass_url_encoding = True,
    )

    data = csv.DictReader(c.fileobj, delimiter = ',')
    result = collections.defaultdict(set)

    for rec in data:

        uniprots = mapping.map_name(
            rec['Gene Symbol'],
            'genesymbol',
            'uniprot',
        )

        for uniprot in uniprots:

            result[uniprot].add(
                CancerGeneCensusAnnotation(
                    tier = int(rec['Tier']),
                    hallmark = rec['Hallmark'].strip().lower() == 'yes',
                    somatic = rec['Somatic'].strip().lower() == 'yes',
                    germline = rec['Germline'].strip().lower() == 'yes',
                    tumour_types_somatic = (
                        multi_field(rec['Tumour Types(Somatic)'])
                    ),
                    tumour_types_germline = (
                        multi_field(rec['Tumour Types(Germline)'])
                    ),
                    cancer_syndrome = (
                        multi_field(rec['Cancer Syndrome'])
                    ),
                    tissue_type = (
                        multi_field(rec['Tissue Type'].replace(' ', ''))
                    ),
                    genetics = rec['Molecular Genetics'].strip() or None,
                    role = (
                        multi_field(rec['Role in Cancer'])
                    ),
                    mutation_type = (
                        multi_field(rec['Mutation Types'])
                    ),
                )
            )

    return dict(result)


def intogen_annotations():
    """
    Returns a list of cancer driver genes according to the IntOGen database.
    """

    IntogenAnnotation = collections.namedtuple(
        'IntogenAnnotation',
        [
            'type',
            'role',
            'curated',
            'oncodrive_role_prob',
        ],
    )


    url = urls.urls['intogen']['drivers_url']

    c = curl.Curl(
        url,
        large = True,
        silent = False,
        files_needed = ['Drivers_type_role.tsv'],
    )

    for _ in xrange(7):

        __ = c.result['Drivers_type_role.tsv'].readline()

    data = csv.DictReader(
        c.result['Drivers_type_role.tsv'],
        delimiter = '\t',
    )
    result = collections.defaultdict(set)

    for rec in data:

        uniprots = mapping.map_name(
            rec['geneHGNCsymbol'],
            'genesymbol',
            'uniprot',
        )

        for uniprot in uniprots:
            
            role_prob, curated = (
                (
                    1.0,
                    True,
                )
                if rec['OncodriveROLE_prob'] == 'Manually curated' else
                (
                    common.float_or_nan(rec['OncodriveROLE_prob']),
                    False,
                )
            )
            
            result[uniprot].add(
                IntogenAnnotation(
                    type = rec['Driver_type'],
                    role = rec['Role'],
                    curated = curated,
                    oncodrive_role_prob = role_prob,
                )
            )

    return result


def get_matrixdb(organism = 9606):

    url = urls.urls['matrixdb']['url']
    c = curl.Curl(url, silent = False, large = True)
    f = c.result
    i = []
    lnum = 0
    for l in f:
        if lnum == 0:
            lnum += 1
            continue

        l = l.replace('\n', '').replace('\r', '')
        l = l.split('\t')
        specA = 0 if l[9] == '-' else int(l[9].split(':')[1].split('(')[0])
        specB = 0 if l[10] == '-' else int(l[10].split(':')[1].split('(')[0])
        if organism is None or (specA == organism and specB == organism):
            pm = [
                p.replace('pubmed:', '') for p in l[8].split('|')
                if p.startswith('pubmed:')
            ]
            met = [
                m.split('(')[1].replace(')', '') for m in l[6].split('|')
                if '(' in m
            ]
            l = [l[0], l[1]]
            interaction = ()
            for ll in l:
                ll = ll.split('|')
                uniprot = ''
                for lll in ll:
                    nm = lll.split(':')
                    if nm[0] == 'uniprotkb' and len(nm[1]) == 6:
                        uniprot = nm[1]
                interaction += (uniprot, )
            interaction += ('|'.join(pm), '|'.join(met))
            if len(interaction[0]) > 5 and len(interaction[1]) > 5:
                i.append(list(interaction))
        lnum += 1
    f.close()
    return i


def get_innatedb(organism = 9606):

    url = urls.urls['innatedb']['url']
    c = curl.Curl(url, silent = False, large = True)
    f = c.result
    i = []
    lnum = 0
    for l in f:
        if lnum == 0:
            lnum += 1
            continue

        l = l.replace('\n', '').replace('\r', '')
        l = l.split('\t')
        specA = 0 if l[9] == '-' else int(l[9].split(':')[1].split('(')[0])
        specB = 0 if l[10] == '-' else int(l[10].split(':')[1].split('(')[0])
        if organism is None or (specA == organism and specB == organism):
            pm = l[8].replace('pubmed:', '')
            l = [l[4], l[5]]
            interaction = ()
            for ll in l:
                ll = ll.split('|')
                hgnc = ''
                uniprot = ''
                for lll in ll:
                    nm = lll.split(':')
                    if nm[0] == 'hgnc':
                        hgnc = nm[1].split('(')[0]
                    if nm[0] == 'uniprotkb' and len(nm[1]) == 6:
                        uniprot = nm[1]
                interaction += (uniprot, hgnc)
            interaction += (pm, )
            i.append(interaction)
        lnum += 1
    f.close()
    s = ''
    for l in i:
        line = ';'.join(list(l)) + "\n"
        if len(line) > 12:
            s += line
    return i


def mitab_field_list(field):
    return common.uniq_list(
        map(lambda x: x.split('(')[1][:-1], field.split('|')))


def mitab_field_uniprot(field):
    uniprots = list(
        filter(lambda x: len(x) == 2 and x[0] == 'uniprotkb',
               map(lambda x: x.split(':'), field.split('|'))))
    if len(uniprots) > 0:
        return uniprots[0][1]
    else:
        return None


def get_dip(url = None,
            organism = 9606,
            core_only = True,
            direct_only = True,
            small_scale_only = True):

    strDipCore = 'dip-quality-status:core'
    strDirect = 'direct interaction'
    strPhysInt = 'physical interaction'
    strSmallS = 'small scale'
    url = urls.urls['dip']['url'] % ('CR' if core_only else '') \
        if url is None else url
    c = curl.Curl(url, silent = False, large = True)
    f = c.result
    i = []
    lnum = 0
    for l in f:
        if lnum == 0:
            lnum += 1
            continue
        l = l.replace('\n', '').replace('\r', '')
        l = l.split('\t')
        specA = int(l[9].split(':')[1].split('(')[0])
        specA = 0 if l[9] == '-' else int(l[9].split(':')[1].split('(')[0])
        specB = 0 if l[10] == '-' else int(l[10].split(':')[1].split('(')[0])
        intProp = mitab_field_list(l[11])
        intProc = mitab_field_list(l[15])
        dipLinkId = l[13]
        expEv = mitab_field_list(l[6])
        conf = l[14]
        if organism is None or (specA == organism and specB == organism):
            if (not core_only or strDipCore in conf) and \
                (not direct_only or strDirect in intProp or
                    strPhysInt in intProp) and \
                    (not small_scale_only or strSmallS in intProc):
                pm = l[8].replace('pubmed:', '').split('|')
                pm = [p for p in pm if not p.startswith('DIP')]
                l = [l[0], l[1]]
                uniprotA = mitab_field_uniprot(l[0])
                uniprotB = mitab_field_uniprot(l[1])
                if uniprotA is not None and uniprotB is not None:
                    i.append([
                        uniprotA, uniprotB, ';'.join(pm), ';'.join(intProp),
                        ';'.join(expEv), dipLinkId
                    ])
        lnum += 1
    f.close()
    return i


def dip_login(user, passwd):
    """
    This does not work for unknown reasons.

    In addition, the binary_data parameter of Curl().__init__() has been changed,
    below updates are necessary.
    """
    bdr = '---------------------------8945224391427558067125853467'
    useragent = 'Mozilla/5.0 (X11; U; Linux i686; en-US; rv:43.0) '\
        'Gecko/20110304 Firefox/43.0'
    loginfname = os.path.join(settings.get('cachedir'), 'dip.logindata.tmp')
    url = urls.urls['dip']['login']
    req_hdrs = ['User-Agent: %s' % useragent]
    c = curl.Curl(
        url,
        cache = False,
        write_cache = False,
        req_headers = req_hdrs,
        return_headers = True,
        debug = True)
    res = c.result
    hdr = c.resp_headers
    cookie = hdr['set-cookie'].split(';')[0]
    cookie2 = '%s%u' % (cookie[:-1], int(cookie[-1]) + 1)
    othercookie = 'DIPID=11133%3A'
    req_hdrs = [
        'Content-type: multipart/form-data; '
        'boundary = %s' % bdr,
        'Cookie: %s' % cookie2,
        'Referer: %s' % url,
        'User-Agent: %s' % useragent,
        'Connection: keep-alive',
    ]
    post = {'lgn': '1', 'login': user, 'pass': passwd, 'Login': 'Login'}
    login = '--%s\r\n\r\nContent-Disposition: form-data; name = "lgn"\r\n\r\n1'\
        '\r\n--%s\r\n\r\nContent-Disposition: form-data; name = "login"\r\n\r\n'\
        '%s\r\n--%s\r\n\r\nContent-Disposition: form-data; name = "pass"\r\n\r'\
        '\n%s\r\n--%s\r\n\r\nContent-Disposition: form-data; name = "Login"\r\n'\
        '\r\nLogin\r\n%s--\r\n' % (bdr, bdr, user, bdr, passwd, bdr, bdr)
    # login = login.replace('\r', '')
    with codecs.open(loginfname, encoding = 'ISO-8859-1', mode = 'w') as f:
        f.write(login)
    c = curl.Curl(
        url,
        cache = False,
        write_cache = False,
        follow = True,
        req_headers = req_hdrs,
        timeout = 10,
        binary_data = loginfname,
        return_headers = True,
        debug = True)
    res = c.result
    hdr = c.resp_headers
    return res, hdr


def spike_interactions(high_confidence = True):
    url = urls.urls['spike']['url']
    c = curl.Curl(
        url, silent = False, large = True, files_needed = ['LatestSpikeDB.xml'])
    spikexml = c.result

    xml = ET.parse(spikexml['LatestSpikeDB.xml'])

    xmlroot = xml.getroot()

    # iterating genes

    bblock = xmlroot.find('BuildingBlock')
    rblock = xmlroot.find('RegulationBlock')
    iblock = xmlroot.find('InteractionBlock')

    genes = {}

    #TODO: Denes, do we need the out var below?
    #
    #out = 'Entrez_A\tGeneSymbol_A\tEntrez_B\tGeneSymbol_B\tIsDirected\tPMID'\
    #    '\tConfidence\tEffect\tAssayType\tDataSource\tDescription\tMechanism\n'

    result = []

    for gene in bblock.findall('Gene'):
        sy = '' if 'name' not in gene.attrib else gene.attrib['name']
        genes[gene.attrib['id']] = (gene.find('XRef').attrib['id'], sy)

    for reg in rblock.findall('Regulation'):
        ds = reg.attrib['dataSource']
        itg = reg.attrib['integrity']
        eff = reg.attrib['effect']
        mec = reg.attrib['mechanism']
        src = reg.find('Source').attrib['ref']
        tgt = reg.find('PhysicalTarget').attrib['ref']
        dcd = "1"
        dsc = '' if reg.find('Description') is None else reg.find('Description').text.\
            replace('\n', ' ')
        asy = '' if 'biologicalAssay' not in reg.attrib else reg.attrib[
            'biologicalAssay']
        refs = reg.findall('Reference')
        pmids = []
        for r in refs:
            pmids.append(r.attrib['pmid'])
        if src in genes and tgt in genes:
            if itg == '1' or not high_confidence:
                result.append([
                    genes[src][0], genes[src][1], genes[tgt][0], genes[tgt][1],
                    dcd, ';'.join(pmids), itg, eff, asy, ds, dsc, mec
                ])

    for ict in iblock.findall('Interaction'):
        ds = ict.attrib['dataSource']
        itg = ict.attrib['integrity']
        eff = '' if 'effect' not in ict.attrib else ict.attrib['effect']
        src = ict.find('ProteinA').attrib['ref']
        tgt = ict.find('ProteinB').attrib['ref']
        dcd = "0"
        mec = ""
        dsc = '' if ict.find('Description') is None else ict.find('Description').text.\
            replace('\n', ' ')
        asy = '' if 'biologicalAssay' not in ict.attrib else ict.attrib[
            'biologicalAssay']
        refs = ict.findall('Reference')
        pmids = []
        for r in refs:
            pmids.append(r.attrib['pmid'])
        if src in genes and tgt in genes:
            if itg == '1' or not high_confidence:
                result.append([
                    genes[src][0], genes[src][1], genes[tgt][0], genes[tgt][1],
                    dcd, ';'.join(pmids), itg, eff, asy, ds, dsc, mec
                ])

    return result


def mppi_interactions(organism = 9606):
    url = urls.urls['mppi']['url']
    c = curl.Curl(url, silent = False, large = True)
    xmlfile = c.gzfile

    prefix = '{net:sf:psidev:mi}'

    result = []

    xml = ET.parse(xmlfile)

    xmlroot = xml.getroot()

    ilist = xmlroot[0][1]

    proteinInteractor = './/%sproteinInteractor' % prefix
    _organism = './/%sorganism' % prefix
    organism = '%u' % organism
    ncbiTaxId = 'ncbiTaxId'
    primaryRef = './/%sprimaryRef' % prefix
    bibref = './/%sbibref' % prefix
    interactionDetection = './/%sinteractionDetection' % prefix
    shortLabel = './/%sshortLabel' % prefix
    fullName = './/%sfullName' % prefix

    #out = 'PMID\tMethods\tUniProt_A\tSwissP/TrEMBL_A\tGene_names_A\t'\
    #'NCBI_tax_ID_A\tUniProt_B\tSwissP/TrEMBL_B\tGene_names_B\tNCBI_tax_ID_B'

    for i in ilist:
        _proteins = i.findall(proteinInteractor)
        if len(_proteins) == 2 and (
                organism is None or
            (_proteins[0].findall(_organism)[0].attrib[ncbiTaxId] == organism
             and _proteins[1].findall(_organism)[0].attrib[ncbiTaxId] ==
             organism)):
            pmids = []
            pms = i.findall(bibref)[0].findall(primaryRef)
            for pm in pms:
                if 'id' in pm.attrib:
                    pmids.append(pm.attrib['id'])
            meths = []
            dets = i.findall(interactionDetection)[0]\
                .findall(shortLabel)
            for m in dets:
                meths.append(m.text)
            proteins = []
            for prot in _proteins:
                thisP = {}
                if 'id' in prot.findall(primaryRef)[0].attrib:
                    thisP['u'] = prot.findall(primaryRef)[0].attrib['id']
                else:
                    thisP['u'] = ''
                thisP['nt'] = prot.findall(primaryRef)[0].attrib['db']
                thisP['gn'] = prot.findall(fullName)[0].text
                thisP['o'] = prot.findall(_organism)[0].attrib[ncbiTaxId]
                proteins.append(thisP)

            result.append([
                ';'.join(pmids), ';'.join(pmids), proteins[0]['u'],
                proteins[0]['nt'], proteins[0]['gn'], proteins[0]['o'],
                proteins[1]['u'], proteins[1]['nt'], proteins[1]['gn'],
                proteins[1]['o']
            ])
    return result


def negatome_pairs():
    url = urls.urls['negatome']['manual']
    c = curl.Curl(url, silent = False, large = True)
    f = c.result
    result = []
    for l in f:
        l = l.strip().split('\t')
        if len(l) == 4:
            l[3] = ';'.join(
                map(lambda x: x.split('-')[1].strip(),
                    filter(lambda x: '-' in x, l[3].replace('–', '-').split(
                        ','))))
        l[0] = l[0].split('-')[0]
        l[1] = l[1].split('-')[0]
        result.append(l)
    return result


def trim_macrophage_gname(gname):
    gname = re.sub(r'\[.*\]', '', re.sub(r'\(.*\)', '', gname))
    gname = re.sub(r'[A-Z]{0,1}[a-z]{1,}', '', gname)
    gname = gname.split(':')
    for i, g in enumerate(gname):
        gname[i] = gname[i].strip()
    return gname


def macrophage_interactions():
    url = urls.urls['macrophage']['url']
    c = curl.Curl(url, silent = False, large = True)
    fname = c.fileobj.name
    del c
    tbl = inputs_common.read_xls(fname)[5:]
    types = ["Protein", "Complex"]
    lst = []
    lnum = 0
    for l in tbl:
        null = ['', '-']
        if len(l) > 11:
            if l[3].strip() in types and l[7].strip() in types:
                alist = trim_macrophage_gname(l[1])
                blist = trim_macrophage_gname(l[5])
                if len(alist) > 0 and len(blist) > 0:
                    for i in alist:
                        for j in blist:
                            if i != j not in null and i not in null:
                                pm = l[11].replace(',',
                                                   '').strip().split('.')[0]
                                if not pm.startswith('INF'):
                                    d = "0" if l[9].strip(
                                    ) == "Binding" else "1"
                                    lst.append([
                                        i, j, l[9].strip(), d, l[10].strip(),
                                        pm
                                    ])
        lnum += 1

    return lst


def get_string_effects(ncbi_tax_id = 9606,
                       stimulation = ['activation'],
                       inhibition = ['inhibition'],
                       exclude = ['expression'],
                       score_threshold = 0):

    effects = []
    if type(stimulation) is list:
        stimulation = set(stimulation)
    if type(inhibition) is list:
        inhibition = set(inhibition)
    if type(exclude) is list:
        exclude = set(exclude)
    url = urls.urls['string']['actions'] % ncbi_tax_id
    c = curl.Curl(url, silent = False, large = True)
    _ = next(c.result)

    for l in c.result:

        l = l.decode('ascii').split('\t')
        if len(l) and l[4] == '1' \
                and int(l[5]) >= score_threshold:
            eff = '+' if l[2] in stimulation \
                else '-' if l[2] in inhibition \
                else '*' if l[2] not in exclude \
                else None
            if eff is not None:
                effects.append([l[0][5:], l[1][5:], eff])

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


def mir2disease_interactions():

    url = urls.urls['mir2dis']['url_rescued']
    c = curl.Curl(url, silent = True, large = True, encoding = 'iso-8859-1')

    return [
        l.strip().split('\t')
        for l in itertools.islice(c.result, 3, None)
    ]

def mirdeathdb_interactions():

    url = urls.urls['mirdeathdb']['url_rescued']
    c = curl.Curl(url, silent = False, large = True)

    _ = next(c.result)

    for l in c.result:

        l = l.strip().split('\t')

        if len(l) < 11:
            continue

        mirnas = l[2].replace('"', '').split(',')
        organism = int(l[9])
        pubmed = l[8]
        geneid = l[10]
        function = '%s_%s' % (l[4], l[5])

        for mirna in mirnas:

            yield (mirna.strip(), geneid, organism, pubmed, function)


def ncrdeathdb_interactions():
    
    
    NcrdeathdbInteraction = collections.namedtuple(
        'NcrdeathdbInteraction',
        [
            'ncrna',
            'protein',
            'ncrna_type',
            'pathway',
            'effect',
            'pmid',
            'organism',
        ]
    )
    
    
    url = urls.urls['ncrdeathdb']['url']
    c = curl.Curl(
        url,
        large = True,
        silent = False,
        encoding = 'iso-8859-1',
    )
    
    data = csv.DictReader(c.fileobj, delimiter = '\t')
    result = []
    
    for rec in data:
        
        typ = rec['RNA Category'].strip()
        
        rna_ids = (
            (rec['miRNA_symbol'],)
                if typ == 'lncRNA' else
            rec['miRBase_ID'].split(',')
        )
        
        for rna_id in rna_ids:
            
            rna_id = rna_id.strip() or None
            protein_id = rec['Gene_Symbol'].strip() or None
            
            if not rna_id and not protein_id:
                
                continue
            
            result.append(NcrdeathdbInteraction(
                ncrna = rna_id,
                protein = protein_id,
                ncrna_type = typ,
                pathway = rec['Pathway'].strip(),
                effect = rec['Action_Mode'].strip() or None,
                pmid = rec['PMID'].strip(),
                organism = int(rec['tax_id'].strip()),
            ))
    
    return result


def mirecords_interactions():

    url = urls.urls['mirecords']['url']
    c = curl.Curl(url, silent = False, large = True)

    tbl = inputs_common.read_xls(c.fileobj.name)

    c.close()

    return (
        (l[6], l[3], l[2], l[1], l[5], l[0].split('.')[0])
        for l in
        ([f.strip() for f in ll] for ll in tbl[1:])
    )


def mirtarbase_interactions():

    url = urls.urls['mirtarbase']['strong']
    c = curl.Curl(url, silent = False, large = True)

    tbl = inputs_common.read_xls(c.fileobj.name)

    c.close()

    for i in xrange(len(tbl)):
        tbl[i][4] = tbl[i][4].split('.')[0]
        tbl[i][8] = tbl[i][8].split('.')[0]

    return tbl[1:]


def lncdisease_interactions():

    url = urls.urls['lncdisease']['url_rescued']
    c = curl.Curl(url, silent = False, large = True)

    for l in c.result:

        l = l.strip().split('\t')

        yield (
            l[1],
            l[2],
            l[3].split('-')[0],
            l[3].split('-')[1] if '-' in l[3] else '',
            l[4].lower(),
            l[6].lower(),
            l[9],
        )


def lncrnadb_interactions():

    renondigit = re.compile(r'[^\d]+')

    url = urls.urls['lncrnadb']['url_rescued']
    c = curl.Curl(
        url,
        silent = False,
        large = True,
        encoding = 'utf-8',
    )

    b = bs4.BeautifulSoup(c.fileobj, 'lxml')

    for res in b.findAll('results'):

        lncrna = res.find('nomenclature').find('name').text

        for sp in res.find('species').findAll('entry'):

            spec = sp.attrs['species'].split('(')[0].strip()

            for assoc in res.find('association').findAll('association'):

                partner  = assoc.find('componentid').text
                typ      = assoc.find('componenttype').text.lower()
                pmid     = renondigit.sub('', assoc.find('pubmedid').text)

                yield (lncrna, partner, typ, spec, pmid)


def transmir_interactions():

    url = urls.urls['transmir']['url']
    c = curl.Curl(
        url,
        silent = False,
        large = True,
        encoding = 'iso-8859-1',
    )
    
    TransmirInteraction = collections.namedtuple(
        'TransmirInteraction',
        [
            'tf_genesymbol',
            'mirna',
            'effect',
            'pubmed',
        ]
    )
    
    result = []

    for l in c.result:

        l = l.strip().split('\t')
        
        result.append(TransmirInteraction(
            tf_genesymbol = l[0],
            mirna = l[1],
            effect = l[4].split('(')[0],
            pubmed = l[5],
        ))
    
    return result


def encode_tf_mirna_interactions():

    url = urls.urls['encode']['tf-mirna']
    c = curl.Curl(url, silent = False, large = True,
                  encoding = 'ascii')

    for l in c.result:

        l = l.strip().split()

        if l[1] == '(TF-miRNA)':

            yield (l[0], l[2])


def _get_imweb():

    def init_fun(resp_hdr):
        return ['Cookie: access-token=%s' % resp_hdr['token']]

    t = int(time.time() * 1000) - 3600000

    loginurl = urls.urls['imweb']['login'] % t

    hdrs = [
        'Host: www.intomics.com',
        'X-Requested-With: XMLHttpRequest',
        'User-Agent: Mozilla/5.0 (X11; Linux x86_64; rv:52.0) Gecko/20100101 Firefox/52.0',
        'Accept-Language: en-US,en;q=0.5',
        'DNT: 1',
        'Connection: keep-alive',
        'Referer: https://www.intomics.com/inbio/map/',
        'Accept: */*'
    ]

    c0 = curl.Curl(loginurl, silent = False, large = False,
                   cache = False, req_headers = hdrs)

    hdrs = hdrs[:-2]

    hdrs.extend([
        'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8',
        'Upgrade-Insecure-Requests: 1',
        'Accept-Encoding: gzip'
    ])

    # 'Host: www.intomics.com' -H 'User-Agent: Mozilla/5.0 (X11; Linux x86_64; rv:52.0) Gecko/20100101 Firefox/52.0' -H 'Accept: text/html,application/xhtml+xml,application/xml;q = 0.9,*/*;q = 0.8' -H 'Accept-Language: en-US,en;q = 0.5' --compressed -H 'Cookie: access_token = '"$token" -H 'DNT: 1' -H 'Connection: keep-alive' -H 'Upgrade-Insecure-Requests: 1'

    hdrs.append('Cookie: access-token=%s' % json.loads(c0.result)['token'])

    url = urls.urls['imweb']['url']

    time.sleep(1)

    c2 = curl.Curl(url, silent = False, large = True,
                   req_headers = hdrs, cache = False,
                   compressed = True)

    return c0, c2


def get_imweb(verbose = 0):

    import pycurl
    import time
    import json

    t = int(time.time() * 1000) - 3600000

    url   = 'https://www.intomics.com/inbio/map/api/'\
            'get_data?file=InBio_Map_core_2016_09_12.tar.gz'
    login = 'https://www.intomics.com/inbio/api/login_guest?ref=&_= %u' % t

    fp_login = open('imweb.login.tmp', 'wb')
    fp_imweb = open('imweb.tmp.tar.gz', 'wb')

    c0 = pycurl.Curl()
    c0.setopt(pycurl.URL, login)
    c0.setopt(pycurl.WRITEFUNCTION, fp_login.write)

    c0.perform()

    fp_login.close()

    with open('imweb.login.tmp', 'r') as fp:
        token = json.loads(fp.read())['token']

    print('Token: %s' % token)

    hdrs = ['Cookie: access-token=%s' % token]

    c1 = pycurl.Curl()
    c1.setopt(pycurl.URL, url)
    c1.setopt(pycurl.WRITEFUNCTION, fp_imweb.write)
    c1.setopt(pycurl.HTTPHEADER, [h.encode('ascii') for h in hdrs])
    c1.setopt(pycurl.VERBOSE, 1)
    c1.setopt(pycurl.DEBUGFUNCTION, print)

    c1.perform()

    fp_imweb.close()


def get_imweb_req():

    import requests
    import time
    import json

    t = int(time.time() * 1000) - 3600000

    url   = 'https://www.intomics.com/inbio/map/api/'\
            'get_data?file=InBio_Map_core_2016_09_12.tar.gz'
    login = 'https://www.intomics.com/inbio/api/login_guest?ref=&_=%u' % t

    r0 = requests.get(login)
    token = json.loads(r0.text)['token']
    hdrs = {'Cookie': 'access-token=%s' % token}

    with open('imweb.tmp.tar.gz', 'wb') as fp:

        r1 = requests.get(url, headers = hdrs, stream = True)

        for block in r1.iter_content(4096):

            fp.write(block)


def get_proteinatlas(normal = True, pathology = True, cancer = True):

    result = {
        'normal':    collections.defaultdict(lambda: {}),
        'pathology': collections.defaultdict(lambda: {})
    }

    def line(l):

        return l.strip('\n\r').split('\t')

    if normal:

        c = curl.Curl(urls.urls['proteinatlas']['normal'],
                    silent = False, large = True)
        fp = list(c.result.values())[0]
        hdr = line(fp.readline().strip())

        for l in fp:

            l = line(l)

            uniprots = mapping.map_name(l[1], 'genesymbol', 'uniprot')
            tissue = '%s:%s' % (l[2], l[3])

            for u in uniprots:
                result['normal'][tissue][u] = (l[4], l[5].strip())

    if cancer or pathology:

        c = curl.Curl(urls.urls['proteinatlas']['pathology'],
                    silent = False, large = True)
        fp = list(c.result.values())[0]
        hdr = line(fp.readline())

        for l in fp:

            l = line(l)
            uniprots = mapping.map_name(l[1], 'genesymbol', 'uniprot')
            tissue   = l[2]

            values = dict(
                (h, float(l[i + 3]) if '.' in l[i + 3] else int(l[i + 3]))
                for i, h in enumerate(hdr[3:])
                if len(l) and len(l[i + 3].strip())
            )

            for u in uniprots:
                result['pathology'][tissue][u] = values

    return result


def proteinatlas_annotations(normal = True, pathology = True, cancer = True):

    LEVELS = ('Not detected', 'Low', 'Medium', 'High')

    ProteinatlasAnnotation = collections.namedtuple(
        'ProtainatlasAnnotation',
        [
            'organ',
            'tissue',
            'level',
            'status',
            'n_not_detected',
            'n_low',
            'n_medium',
            'n_high',
            'prognostic',
            'favourable',
            'score',
            'pathology',
        ],
    )
    ProteinatlasAnnotation.__new__.__defaults__ = (
        (None,) * 4 + (False, False, None, False)
    )


    def n_or_none(ex, key):

        return ex[key] if key in ex else None


    data = get_proteinatlas(
        normal = normal,
        pathology = pathology,
        cancer = cancer,
    )

    result = collections.defaultdict(set)

    if normal:

        for tissue, gex in iteritems(data['normal']):
            
            organ = tissue
            
            if ':' in tissue:
                
                organ, tissue = tissue.split(':')
            
            organ = organ.strip()
            tissue = tissue.strip()

            for uniprot, ex in iteritems(gex):

                uniprots = mapping.map_name(uniprot, 'uniprot', 'uniprot')
                
                for _uniprot in uniprots:
                    
                    result[_uniprot].add(
                        ProteinatlasAnnotation(
                            organ = organ,
                            tissue = tissue,
                            level = ex[0],
                            status = ex[1],
                        )
                    )

    if pathology or cancer:

        for condition, gex in iteritems(data['pathology']):

            for uniprot, ex in iteritems(gex):

                try:
                    effect, score = next(
                        i for i in iteritems(ex) if i[0] not in LEVELS
                    )
                    prognostic = not effect.startswith('unprognostic')
                    favourable = not effect.endswith('unfavourable')

                except StopIteration:

                    prognostic, favourable, score = None, None, None
                
                uniprots = mapping.map_name(uniprot, 'uniprot', 'uniprot')
                
                for _uniprot in uniprots:
                    
                    result[_uniprot].add(
                        ProteinatlasAnnotation(
                            organ = condition,
                            tissue = condition,
                            level = max(
                                (i for i in iteritems(ex) if i[0] in LEVELS),
                                key = lambda i: i[1],
                                default = (None,),
                            )[0],
                            status = None,
                            n_not_detected = n_or_none(ex, 'Not detected'),
                            n_low = n_or_none(ex, 'Low'),
                            n_medium = n_or_none(ex, 'Medium'),
                            n_high = n_or_none(ex, 'High'),
                            prognostic = prognostic,
                            favourable = favourable,
                            score = score,
                            pathology = True,
                        )
                    )

    return result


def proteinatlas_subcellular_annotations():
    
    
    ProteinatlasSubcellularAnnotation = collections.namedtuple(
        'ProteinatlasSubcellularAnnotation',
        [
            'location',
            'status',
        ],
    )
    
    
    url = urls.urls['proteinatlas']['subcell']
    
    c = curl.Curl(
        url,
        large = True,
        silent = False,
        default_mode = 'r',
    )
    reader = csv.DictReader(
        c.files_multipart['subcellular_location.tsv'],
        delimiter = '\t',
    )
    
    result = collections.defaultdict(set)
    
    for rec in reader:
        
        uniprots = mapping.map_name(rec['Gene name'], 'genesymbol', 'uniprot')
        
        for uniprot in uniprots:
            
            for status in ('Enhanced', 'Supported', 'Uncertain'):
                
                if not rec[status]:
                    
                    continue
                
                for location in rec[status].split(';'):
                    
                    result[uniprot].add(ProteinatlasSubcellularAnnotation(
                        location = location,
                        status = status,
                    ))
    
    return dict(result)


def proteinatlas_secretome_annotations():
    
    
    ProteinatlasSecretomeAnnotation = collections.namedtuple(
        'ProteinatlasSecretomeAnnotation',
        [
            'mainclass',
            'secreted',
        ],
    )
    
    
    url = urls.urls['proteinatlas']['secretome']
    
    c = curl.Curl(
        url,
        large = True,
        silent = False,
    )
    
    reader = inputs_common.read_xls(c.fileobj.name)[1:]
    result = collections.defaultdict(set)
    
    for rec in reader:
        
        for uniprot_original in rec[2].split(','):
            
            uniprots = mapping.map_name(
                uniprot_original,
                'uniprot',
                'uniprot',
            )
            
            for uniprot in uniprots:
                
                result[uniprot].add(ProteinatlasSecretomeAnnotation(
                    mainclass = rec[3],
                    secreted = 'secreted' in rec[3].lower(),
                ))
    
    return dict(result)


def get_tfregulons_old(
        levels = {'A', 'B'},
        only_curated = False
    ):
    """
    Retrieves TF-target interactions from TF regulons.

    :param set levels:
        Confidence levels to be used.
    :param bool only_curated:
            Retrieve only literature curated interactions.

    Details
    -------
    TF regulons is a comprehensive resource of TF-target interactions
    combining multiple lines of evidences: literature curated databases,
    ChIP-Seq data, PWM based prediction using HOCOMOCO and JASPAR matrices
    and prediction from GTEx expression data by ARACNe.

    For details see https://github.com/saezlab/DoRothEA.
    """

    url = urls.urls['tfregulons']['url'] % (
        'all' if 'E' in levels else
        'ABCD' if 'D' in levels else
        'ABC' if 'C' in levels else
        'AB' if 'B' in levels else
        'A'
    )

    c = curl.Curl(url, silent = False, large = True)
    _ = next(c.result)

    return (
        list(
            itertools.chain(
                ll[:4],
                (s == 'TRUE' for s in ll[4:8]),
                ll[-4:],
                [','.join(s for s in ll[-4:] if s)]
                if not only_curated else ll[8]
            )
        )
        for ll in (
            l.strip('\n\r').split('\t') for l in c.result
        ) if (
            ll[3] in levels and
            not only_curated or ll[4] == 'TRUE'
        )
    )


def get_tfregulons(
        levels = {'A', 'B'},
        only_curated = False
    ):
    """
    Retrieves TF-target interactions from TF regulons.

    :param set levels:
        Confidence levels to be used.
    :param bool only_curated:
            Retrieve only literature curated interactions.

    Details
    -------
    TF regulons is a comprehensive resource of TF-target interactions
    combining multiple lines of evidences: literature curated databases,
    ChIP-Seq data, PWM based prediction using HOCOMOCO and JASPAR matrices
    and prediction from GTEx expression data by ARACNe.

    For details see https://github.com/saezlab/DoRothEA.
    """
    
    
    def process_sources(sources):
        
        upper = (
            'jaspar',
            'trred',
            'kegg',
            'trrust',
            'tred',
            'trrd',
            'hocomoco',
            'fantom4',
        )
        
        special_case = {
            'tfact': 'TFactS',
            'oreganno': 'ORegAnno',
            'reviews': 'DoRothEA-reviews',
            'HOCOMOCO_v11': 'HOCOMOCO-v11',
            'JASPAR_v2018': 'JASPAR-v2018',
        }
        
        
        revia = re.compile(r',|_via_')
        
        sources = functools.reduce(
            lambda s, r: s.replace(r, r.upper()),
            upper,
            sources,
        )
        
        sources = functools.reduce(
            lambda s, r: s.replace(*r),
            iteritems(special_case),
            sources,
        )
        
        return ','.join(
            '%s_DoRothEA' % s
            for s in revia.split(sources)
        )
    
    
    evidence_types = (
        'chipSeq',
        'TFbindingMotif',
        'coexpression',
        'curateddatabase'
    )
    
    DorotheaInteraction = collections.namedtuple(
        'DorotheaInteraction',
        [
            'tf',
            'target',
            'effect',
            'level',
            'curated',
            'chipseq',
            'predicted',
            'coexp',
            'curated_sources',
            'chipseq_sources',
            'predicted_sources',
            'coexp_sources',
            'all_sources',
            'pubmed',
            'kegg_pathways',
        ]
    )
    
    url = urls.urls['tfregulons_git']['url']

    c = curl.Curl(
        url,
        silent = False,
        large = True,
        files_needed = ['database.csv'],
    )

    reader = csv.DictReader(c.result['database.csv'])

    for rec in reader:

        # process only the ones of the requested levels or if curated
        if (
            rec['score'] not in levels and
            not (
                only_curated and
                rec['is_evidence_curateddatabase'] == 'TRUE'
            )
        ):

            continue
        
        rec = dict(
            (k, v if v not in {'-', 'none'} else '')
            for k, v in iteritems(rec)
        )
        
        yield DorotheaInteraction(
            **dict(zip(
                DorotheaInteraction._fields,
                itertools.chain(
                    # TF, target, effect, score
                    (
                        rec[key] for key in 
                        ('TF', 'target', 'effect', 'score')
                    ),
                    # boolean values for curated, chipseq, motif pred. and coexp
                    (
                        rec['is_evidence_%s' % key] == 'TRUE'
                        for key in evidence_types
                    ),
                    # databases & datasets
                    (
                        rec['which_%s' % key]
                        for key in evidence_types
                    ),
                    # all data sources (or only the curated ones)
                    (
                        process_sources(
                            ','.join(
                            rec[key]
                                for key in
                                ('which_%s' % evt for evt in evidence_types)
                                if rec[key]
                            )
                                if not only_curated else
                            rec['which_curateddatabase']
                        ),
                    ),
                    # PubMed and KEGG pw
                    (
                        rec['pubmedID_from_curated_resources'],
                        rec['kegg_pathway'],
                    )
                )
            ))
        )


# synonyms
tfregulons_interactions = get_tfregulons
get_dorothea = get_tfregulons
dorothea_interactions = get_tfregulons


def stitch_interactions(threshold = None):

    url = urls.urls['stitch']['actions']

    c = curl.Curl(url, silent = False, large = True)

    _ = next(c.result)

    sep = re.compile(r'[sm\.]')

    for l in c.result:

        l = l.decode('utf-8').strip().split('\t')

        score = int(l[5])

        if threshold is not None and score < threshold:

            continue

        try:

            a = sep.split(l[0])[1]
            b = sep.split(l[1])[1]

        except IndexError:

            print(l[1])

        if l[4] == 'f':

            a, b = b, a

        yield a, b, l[2], l[3], int(l[5])

def get_cspa(organism = 9606):

    sheets = {
        'Human': 'Table A',
        'Mouse': 'Table B',
    }

    str_organism = taxonomy.taxids[organism].capitalize()

    url = urls.urls['cspa']['url']
    c = curl.Curl(url, large = True, silent = False)
    xlsname = c.fname
    del(c)
    raw = inputs_common.read_xls(xlsname, sheets[str_organism])[1:]

    return mapping.map_names((r[1] for r in raw), 'uniprot', 'uniprot')


def surfaceome_annotations():
    """
    Downloads the "In silico human surfaceome".
    Dict with UniProt IDs as key and tuples of surface prediction score,
    class and subclass as values (columns B, N, S and T of table S3).
    """

    url = urls.urls['surfaceome']['url']
    c = curl.Curl(url, large = True, silent = False)
    xlsname = c.fname
    del(c)
    raw = inputs_common.read_xls(xlsname, 'in silico surfaceome only')[2:]

    return dict(
        (
            uniprot, # uniprot
            (
                float(r[13]), # score
                r[18] if r[18] else None, # class
                set(r[19].split(';')) if r[19] else set(), # subclass
            )
        )
        for r in raw
        for uniprot in mapping.map_name(r[1], 'uniprot', 'uniprot')
    )

def matrisome_annotations(organism = 9606):
    """
    Downloads MatrisomeDB 2.0, a database of extracellular matrix proteins.
    Returns dict where keys are UniProt IDs and values are tuples of
    classes, subclasses and notes.
    """
    
    MatrisomeAnnotation = collections.namedtuple(
        'MatrisomeAnnotation',
        ['mainclass', 'subclass', 'subsubclass']
    )
    
    tax_names = {
        10090: ('Murine', 'mm'),
        9606:  ('Human',  'hs'),
    }

    url = urls.urls['matrisome']['url_xls'] % tax_names[organism]
    c = curl.Curl(url, large = True, silent = False)
    xlsname = c.fname
    del(c)
    raw = inputs_common.read_xls(xlsname)[1:]

    result = collections.defaultdict(set)
    
    for r in raw:
        
        uniprots = set(r[7].split(':'))
        uniprots.discard('')
        
        if not uniprots:
            
            continue
        
        uniprots = mapping.map_names(uniprots, 'uniprot', 'uniprot')
        
        for uniprot in uniprots:
            
            result[uniprot].add(
                MatrisomeAnnotation(
                    mainclass = r[0].strip(),
                    subclass = r[1].strip(),
                    subsubclass = r[10].strip() or None,
                )
            )
    
    return dict(result)


def __matrisome_annotations_2():
    """
    This I made only to find out why certain proteins are missing from this
    output. I will contact Matrisome people to ask why.
    """

    url = urls.urls['matrisome']['url_dl']
    c = curl.Curl(url, large = True, silent = False)

    _ = next(c.result)

    return set(r.split(',')[1] for r in c.result)

def membranome_annotations():

    membr_url = urls.urls['membranome']['baseurl'] % ('membranes', '')
    c = curl.Curl(membr_url, large = True, silent = False)
    membr_data = json.loads(c.fileobj.read())
    del c

    membr = dict((m['id'], m) for m in membr_data['objects'])

    page = 1
    prot_all = []

    prg = progress.Progress(7, 'Downloading Membranome', 1)

    while True:

        prg.step()

        prot_url = urls.urls['membranome']['baseurl'] % (
            'proteins',
            '?pageSize=1000&pageNum=%u' % page,
        )
        c = curl.Curl(prot_url, large = True, silent = True)
        prot = json.loads(c.fileobj.read())

        prot_all.extend(prot['objects'])

        if prot['page_end'] >= prot['total_objects']:

            break

        page = prot['page_num'] + 1

    prg.terminate()

    for p in prot_all:
        
        uniprots = mapping.map_name(p['uniprotcode'], 'uniprot', 'uniprot')
        
        for uniprot in uniprots:
            
            yield (
                uniprot,
                membr[p['membrane_id']]['name'],
                membr[p['membrane_id']]['topology_in']
                    if p['topology_show_in'] else
                membr[p['membrane_id']]['topology_out'],
            )


def opm_annotations(organism = 9606):


    reparentheses = re.compile(r'\((.*)\)')
    regenesymbol  = re.compile(r' ([A-Z0-9]{3,}) ')

    def get_dict(name):

        result = {}
        url = urls.urls['opm'][name]
        c = curl.Curl(url, large = True, silent = False)
        data = csv.DictReader(c.result, delimiter = ',')

        for rec in data:
            result[rec['id']] = rec['name']

        return result


    OpmAnnotation = collections.namedtuple(
        'OpmAnnotation',
        ['membrane', 'family', 'transmembrane'],
    )

    result = collections.defaultdict(set)

    organism_name = (
        taxonomy.phosphoelm_taxids[organism]
            if organism in taxonomy.phosphoelm_taxids else
        None
    )

    types = get_dict('types')
    families = get_dict('families')

    url = urls.urls['opm']['proteins']
    c = curl.Curl(url, silent = False, large = True)

    data = csv.DictReader(c.result, delimiter = ',')

    for rec in data:

        if organism_name and rec['species_name_cache'] != organism_name:

            continue

        name = rec['name']

        names = [
            name,
            name.split('(')[0],
            name.split(',')[0],
        ]

        m = reparentheses.search(name)
        if m:
            names.append(m.groups()[0])

        genesymbols = regenesymbol.findall(name)

        for this_name in names:

            uniprot = mapping.map_name0(this_name, 'protein-name', 'uniprot')

            if uniprot:

                break

        if not uniprot:

            for gs in genesymbols:

                uniprot = (
                    mapping.map_name0(this_name, 'genesymbol', 'uniprot')
                )

                if uniprot:

                    break

        if not uniprot:

            continue

        result[uniprot].add(
            OpmAnnotation(
                membrane = rec['membrane_name_cache'],
                family = rec['family_name_cache'],
                transmembrane = types[rec['type_id']] == 'Transmembrane',
            )
        )

    return result


def topdb_annotations(ncbi_tax_id = 9606):

    TopdbAnnotation = collections.namedtuple(
        'TopdbAnnotation',
        ['membrane', 'topology', 'score', 'tmregions'],
    )

    result = collections.defaultdict(set)

    url = urls.urls['topdb']['url']
    c = curl.Curl(
        url,
        large = True,
        default_mode = 'rb',
        silent = False,
    )

    parser = etree.iterparse(c.fileobj, events = ('start', 'end'))

    result = collections.defaultdict(set)
    root = next(parser)
    used_elements = []

    for ev, elem in parser:

        if ev == 'end' and elem.tag == 'TOPDB':

            used_elements.append(elem)

            organism = elem.find('Organism').text

            if (
                organism not in taxonomy.latin_name_to_ncbi_tax_id or
                taxonomy.latin_name_to_ncbi_tax_id[organism] != ncbi_tax_id
            ):

                continue

            tag_uniprots = elem.find('./CrossRef/UniProt')

            if tag_uniprots is None:

                continue

            uniprots = [u.text for u in tag_uniprots.findall('AC')]
            uniprots = set(
                mapping.map_name0(
                    u,
                    'uniprot',
                    'uniprot',
                    ncbi_tax_id = ncbi_tax_id,
                )
                for u in uniprots
            )

            if not uniprots:

                continue

            membranes = set(
                mem
                for tag_mem in elem.findall('Membrane')
                for mem in tag_mem.text.split(';')
            )

            ntm = 0
            score = 0
            topologies = ()
            tag_topo = elem.find('Topology')

            if tag_topo is not None:

                ntm = int(tag_topo.find('Numtm').attrib['Count'])
                score = int(tag_topo.find('Reliability').text)

                topologies = set(
                    tag_reg.attrib['Loc']
                    for tag_reg in tag_topo.findall('./Regions/Region')
                )

            if not membranes:

                membranes = (None,)

            if not topologies:

                topologies = (None,)

            for topology, membrane, uniprot in itertools.product(
                topologies,
                membranes,
                uniprots,
            ):

                result[uniprot].add(
                    TopdbAnnotation(
                        membrane = membrane,
                        topology = topology,
                        tmregions = ntm,
                        score = score,
                    )
                )

        # removing used elements to keep memory low
        if len(used_elements) > 2000:

            for _ in xrange(1000):

                e = used_elements.pop(0)
                e.clear()

    # closing the XML
    c.fileobj.close()
    del c

    return result


def hgnc_genegroups():

    HGNCGeneGroupAnnotation = collections.namedtuple(
        'HGNCGeneGroupAnnotation',
        ['mainclass'],
    )
    result = collections.defaultdict(set)

    url = urls.urls['hgnc']['groups']
    c = curl.Curl(url, large = True, silent = False)

    _ = next(c.result)

    for rec in c.result:

        rec = rec.split('\t')
        uniprots = {u.strip() for u in rec[2].split(',')}
        uniprots.discard('')
        
        if not uniprots:
            
            continue
        
        uniprots = mapping.map_names(uniprots, 'uniprot', 'uniprot')

        if not uniprots:

            continue

        groups = rec[3].split('|')

        for group in groups:

            group = group.strip()

            if group:
                
                for uniprot in uniprots:
                    
                    result[uniprot].add(
                        HGNCGeneGroupAnnotation(mainclass = group)
                    )

    return result


def zhong2015_annotations():
    """
    From 10.1111/nyas.12776 (PMID 25988664).
    """

    types = {
        'i': 'iCAM',
        'm': 'matrix adhesion',
        'ag': 'axonal guidance',
        'aj': 'adherens junction',
        'c': 'cell-cell adhesion',
        'fa': 'focal adhesion',
        'tj': 'tight junction',
        'my': 'myelin interactions',
    }

    Zhong2015Annotation = collections.namedtuple(
        'Zhong2015Annotation',
        ['type'],
    )
    result = collections.defaultdict(set)

    url = urls.urls['zhong2015']['url']
    c = curl.Curl(url, silent = False, large = True)

    _ = next(c.result)

    for rec in c.result:

        rec = rec.strip().split('\t')

        uniprot = mapping.map_name0(rec[0], 'genesymbol', 'uniprot')

        if uniprot:

            result[uniprot].add(
                Zhong2015Annotation(type = types[rec[2]])
            )

    return result


def adhesome_annotations():

    AdhesomeAnnotation = collections.namedtuple(
        'AdhesomeAnnotation',
        ['mainclass', 'intrinsic'],
    )

    result = collections.defaultdict(set)

    url = urls.urls['adhesome']['components']
    c = curl.Curl(url, large = True, silent = False)

    data = csv.DictReader(c.result, delimiter = ',')

    for rec in data:

        uniprots = rec['Swiss-Prot ID']

        for uniprot in uniprots.split(','):

            uniprot = uniprot.strip()

            if uniprot == 'null':

                continue
            
            for _uniprot in mapping.map_name(uniprot, 'uniprot', 'uniprot'):

                result[uniprot].add(AdhesomeAnnotation(
                    mainclass = rec['Functional Category'].strip(),
                    intrinsic = rec['FA'].strip() == 'Intrinsic Proteins',
                ))

    return result


def get_exocarta(organism = 9606, types = None):
    """
    :param set types:
        Molecule types to retrieve. Possible values: `protein`, `mrna`.
    """

    return _get_exocarta_vesiclepedia(
        database = 'exocarta',
        organism = organism,
        types = types,
    )

def get_vesiclepedia(organism = 9606, types = None):
    """
    :param set types:
        Molecule types to retrieve. Possible values: `protein`, `mrna`.
    """

    return _get_exocarta_vesiclepedia(
        database = 'vesiclepedia',
        organism = organism,
        types = types,
    )

def _get_exocarta_vesiclepedia(
        database = 'exocarta',
        organism = 9606,
        types = None
    ):
    """
    :param str database:
        Which database to download: ExoCarta or Vesiclepedia.
    :param set types:
        Molecule types to retrieve. Possible values: `protein`, `mrna`.
    """

    database = database.lower()

    types = types or {'protein'}

    organism = taxonomy.phosphoelm_taxids[organism]

    taxid_rev = dict((v, k) for k, v in iteritems(taxonomy.phosphoelm_taxids))

    # collecting the references
    url_s = urls.urls[database]['url_study']
    c = curl.Curl(url_s, large = True, silent = False)
    _ = next(c.result)

    studies = {}

    for s in c.result:

        s = s.split('\t')

        organisms = tuple(
            taxid_rev[t.strip()]
            for t in s[2].split('|')
            if t.strip() in taxid_rev
        )

        if not organisms:

            continue

        stud = (
            s[1] if s[1] != '0' else None, # PubMed ID
            organisms, # organism
            s[4], # sample source (cell type, tissue)
        )

        if database == 'vesiclepedia':

            vtype = s[11].strip()

            stud += (
                tuple(vtype.split('/')) if vtype else (),
            )

        studies[int(s[0])] = tuple(stud)

    # processing proteins
    url_p = urls.urls[database]['url_protein']
    c = curl.Curl(url_p, large = True, silent = False)
    _ = next(c.result)

    for s in c.result:

        s = s.split('\t')

        if s[4] != organism or s[1] not in types:

            continue

        yield (
            s[2], # Entrez ID
            s[3], # Gene Symbol
            taxid_rev[s[4]], # NCBI Taxonomy ID
            studies[int(s[5])], # study reference
        )


def _matrixdb_protein_list(category, organism = 9606):
    """
    Returns a set of proteins annotated by MatrixDB.

    :arg str category:
        The protein annotation category. Possible values: `ecm`, `membrane`
        or `secreted`.
    """

    url = urls.urls['matrixdb']['%s_proteins' % category]
    c = curl.Curl(url, silent = False, large = True)

    proteins = set()

    # header row
    _ = next(c.result)

    for l in c.result:

        if not l:

            continue

        proteins.add(
            l.strip().replace('"', '').split('\t')[0]
        )
    
    proteins = mapping.map_names(proteins, 'uniprot', 'uniprot')

    if organism:

        uniprots = all_uniprots(organism = organism, swissprot = True)
        proteins = proteins & set(uniprots)

    return proteins


def matrixdb_membrane_proteins(organism = 9606):
    """
    Returns a set of membrane protein UniProt IDs retrieved from MatrixDB.
    """

    return _matrixdb_protein_list('membrane', organism = organism)


def matrixdb_secreted_proteins(organism = 9606):
    """
    Returns a set of secreted protein UniProt IDs retrieved from MatrixDB.
    """

    return _matrixdb_protein_list('secreted', organism = organism)


def matrixdb_ecm_proteins(organism = 9606):
    """
    Returns a set of ECM (extracellular matrix) protein UniProt IDs
    retrieved from MatrixDB.
    """

    return _matrixdb_protein_list('ecm', organism = organism)


def matrixdb_annotations(organism = 9606):

    MatrixdbAnnotation = collections.namedtuple(
        'MatrixdbAnnotation',
        ('mainclass',),
    )
    annot = collections.defaultdict(set)

    for cls in ('membrane', 'secreted', 'ecm'):

        cls_annot = MatrixdbAnnotation(mainclass = cls)

        method = globals()['matrixdb_%s_proteins' % cls]

        for uniprot in method(organism = organism):

            annot[uniprot].add(cls_annot)

    return dict(annot)


def get_compartments_localizations(
        organism = 9606,
        literature = True,
        high_throughput = False,
        text_mining = False,
        predictions = False,
    ):

        pass


def locate_localizations(
        organism = 9606,
        literature = True,
        external = True,
        predictions = False,
    ):

    record = collections.namedtuple(
        'LocateAnnotation',
        ('source', 'location', 'cls', 'pmid', 'score'),
    )
    record.__new__.__defaults__ = (None, None, None)

    organism_uniprots = set(
        all_uniprots(organism = organism, swissprot = True)
    )

    organism_str = taxonomy.taxids[organism]
    url = urls.urls['locate']['url'] % organism_str
    fname = url.split('/')[-1][:-4]

    c = curl.Curl(
        url,
        large = True,
        default_mode = 'rb',
        silent = False,
        files_needed = [fname],
    )
    c.result[fname]

    parser = etree.iterparse(c.result[fname], events = ('start', 'end'))

    result = collections.defaultdict(set)
    root = next(parser)
    used_elements = []

    for ev, elem in parser:

        if ev == 'end' and elem.tag == 'LOCATE_protein':

            tag_protein = elem.find('protein')
            this_uniprot = None
            this_entrez  = None
            this_organism = (
                tag_protein.find('organism').text
                    if tag_protein is not None else
                None
            )
            this_class = (
                tag_protein.find('class').text
                    if tag_protein is not None else
                None
            )

            xrefs = elem.find('xrefs')

            if xrefs is None:

                continue

            for xref in xrefs.findall('xref'):

                src = xref.find('source')
                src_name = src.find('source_name').text

                if src_name == 'UniProtKB-SwissProt':

                    this_uniprot = src.find('accn').text

                if src_name == 'Entrez Gene':

                    this_entrez = src.find('accn').text

                if src_name == 'UniProt/SPTrEMBL' and this_uniprot is None:

                    this_uniprot = src.find('accn').text

            # if we don't know what it is, does not make sense to proceed
            if this_uniprot is None and this_entrez is None:

                continue

            if this_uniprot:

                this_uniprots = mapping.map_name(
                    this_uniprot,
                    'uniprot',
                    'uniprot',
                    ncbi_tax_id = organism,
                )

            if not this_uniprots and this_entrez:

                this_uniprots = mapping.map_name(
                    this_entrez,
                    'entrez',
                    'uniprot',
                    ncbi_tax_id = organism,
                )

            this_uniprots = set(this_uniprots) & organism_uniprots

            # if we don't know what it is, does not make sense to proceed
            if not this_uniprots:

                continue

            if external:

                # External database annotations
                extannot = elem.find('externalannot')

                if extannot is not None:

                    for extannotref in extannot.findall('reference'):

                        sources = []

                        for src in extannotref.findall('source'):

                            src_name = src.find('source_name')

                            if src_name is not None:

                                sources.append(src_name.text)

                        sources = ';'.join(sources) if sources else None

                        locations =  extannotref.find('locations')

                        if locations is not None:

                            for location in locations.findall('location'):

                                for loc in location.iterchildren():

                                    if loc.tag[:4] == 'tier':

                                        this_loc = loc.text.lower().split(',')

                                        for uniprot in this_uniprots:

                                            for _loc in this_loc:

                                                result[uniprot].add(record(
                                                    source = sources,
                                                    location = _loc.strip(),
                                                    cls = this_class,
                                                    score = None,
                                                ))

            if predictions:

                # Predictions
                sclpred = elem.find('scl_prediction')

                if sclpred is not None:

                    for sclpred_src in sclpred.findall('source'):

                        score = float(sclpred_src.find('evaluation').text)

                        if score == 0.0:

                            continue

                        this_src = sclpred_src.find('method').text
                        this_loc = sclpred_src.find('location').text.lower()

                        if this_loc == 'no prediction':

                            continue

                        for uniprot in this_uniprots:

                            result[uniprot].add(record(
                                source = this_src,
                                location = this_loc,
                                cls = this_class,
                                score = score,
                            ))

            if literature:

                # Literature curation
                lit = elem.find('literature')

                if lit is not None:

                    for litref in lit.findall('reference'):

                        locs = set()

                        for lloc in litref.find('locations').findall('location'):

                            for loc in lloc.iterchildren():

                                if loc.tag[:4] == 'tier':

                                    locs.add(loc.text.lower())

                        pmid = litref.find('source')
                        pmid = None if pmid is None else pmid.find('accn').text

                        for loc in locs:

                            for uniprot in this_uniprots:

                                result[uniprot].add(record(
                                    source = 'literature',
                                    location = loc,
                                    pmid = pmid,
                                    cls = this_class,
                                    score = None,
                                ))

        used_elements.append(elem)

        # removing used elements to keep memory low
        if len(used_elements) > 1000:

            for _ in xrange(500):

                e = used_elements.pop(0)
                e.clear()

    # closing the XML
    c.fileobj.close()
    del c

    return result


def cancersea_annotations():
    """
    Retrieves genes annotated with cancer funcitonal states from the
    CancerSEA database.
    """
    
    
    CancerseaAnnotation = collections.namedtuple(
        'CancerseaAnnotation',
        [
            'state',
        ],
    )
    
    
    annotations = collections.defaultdict(set)
    
    url = urls.urls['cancersea']['url']
    c = curl.Curl(url, silent = False, large = False)
    soup = bs4.BeautifulSoup(c.result, 'html.parser')
    
    for row in soup.find_all('tbody')[1].find_all('tr'):
        
        state = row.find_all('td')[0].text
        url_end = row.find_all('td')[-1].find('a').attrs['href']
        data_url = urls.urls['cancersea']['data_url'] % url_end
        c = curl.Curl(data_url, silent = False, large = True)
        
        _ = next(c.result)
        
        for line in c.result:
            
            line = line.strip().split('\t')
            
            uniprots = mapping.map_name(line[1], 'genesymbol', 'uniprot')
            
            for uniprot in uniprots:
                
                annotations[uniprot].add(
                    CancerseaAnnotation(state = state)
                )
    
    return dict(annotations)


def get_protmapper():
    """
    Returns the raw records as read by ``csv.DictReader``.
    From Bachman et al. 2019 "Assembling a phosphoproteomic knowledge base
    using ProtMapper to normalize phosphosite information from databases and
    text mining",
    https://www.biorxiv.org/content/10.1101/822668v3.supplementary-material
    """
    
    url = urls.urls['protmapper']['url']
    files = urls.urls['protmapper']['files']
    c = curl.Curl(url, large = True, silent = False, files_needed = files)
    
    evidences = collections.defaultdict(list)
    
    for rec in csv.DictReader(c.files_multipart['evidences.csv']):
        
        evidences[rec['ID']].append(rec)
    
    records = list(csv.DictReader(c.files_multipart['export.csv']))
    
    return records, evidences


def protmapper_ptms(
        only_evidences = None,
        only_literature = False,
        interactions = False,
    ):
    """
    :arg str,set,NoneType only_evidences:
        Keep only the interactions with these evidence type, e.g. `VALID`.
        See the 'descriptions' column in the 'evidences.csv' supplementary
        table.
    """
    
    databases = {
        'signor': 'Signor',
        'psp': 'PhosphoSite',
        'sparser': 'Sparser',
        'reach': 'REACH',
        'pid': 'NCI-PID',
        'reactome': 'Reactome',
        'rlimsp': 'RLIMS-P',
        'bel': 'BEL-Large-Corpus',
    }
    
    result = []
    only_evidences = common.to_set(only_evidences)
    
    records, evidences = get_protmapper()
    
    for rec in records:
        
        if rec['CTRL_NS'] != 'UP':
            
            continue
        
        if only_evidences:
            
            ev_types = {
                ev['DESCRIPTION']
                for ev in evidences[rec['ID']]
            }
            
            if not only_evidences & ev_types:
                
                continue
        
        references = {
            ev['PMID']
            for ev in evidences[rec['ID']]
            if ev['PMID']
        }
        
        if only_literature and not references:
            
            continue
        
        typ = (
            'phosphorylation'
                if rec['CTRL_IS_KINASE'] == 'True' else
            'unknown'
        )
        sources = {
            '%s_ProtMapper' % (
                databases[source] if source in databases else source
            )
            for source in rec['SOURCES'].strip('"').split(',')
        }
        
        if interactions:
            
            result.append([
                rec['CTRL_ID'],
                rec['TARGET_UP_ID'],
                sources,
                references,
            ])
            
        else:
            
            result.append({
                'kinase': rec['CTRL_ID'],
                'resaa': rec['TARGET_RES'],
                'resnum': int(rec['TARGET_POS']),
                'references': references,
                'substrate': rec['TARGET_UP_ID'],
                'databases': sources,
            })
    
    return result


def protmapper_interactions(**kwargs):
    
    return protmapper_ptms(interactions = True, **kwargs)


def msigdb_download(
        registered_email = None,
        collection = 'msigdb',
        id_type = 'symbols',
        force_download = False,
    ):
    """
    Downloads and preprocesses a collection of gmt format gene sets from
    MSigDB. Returns dict of sets with gene set names as keys and molecular
    identifiers as values.
    
    :arg str,NoneType registered_email:
        An email address registered at MSigDB. If `None` the `msigdb_email`
        from ``pypath.settings`` will be used.
    :arg str collection:
        The name of the gene set collection. For available collections (e.g.
        `h.all` or `c2.cpg`) refer to the MSigDB website:
        http://software.broadinstitute.org/gsea/downloads.jsp#msigdb
        The default value `msigdb` contains all the genesets however you
        won't be able to distinguish which geneset comes from which
        collection. For this you need to download the collections one by one.
    :arg str id_type:
        MSigDB provides Gene Symbols (`symbols`) and Entrez Gene IDs
        (`entrez`).
    :arg bool force_download:
        Download even if cache content is available.
    """
    
    registered_email = registered_email or settings.get('msigdb_email')
    
    if not registered_email:
        
        _log(
            'To download MSigDB you must provide an email address '
            'you have previously registered at '
            '`http://software.broadinstitute.org/gsea/register.jsp`. '
            'Could not proceed, returning empty dict.'
        )
        
        return {}
    
    url = urls.urls['msigdb']['url'] % (
        collection,
        id_type,
    )
    
    req_headers_1 = []
    
    c_nocall = curl.Curl(
        url,
        call = False,
        process = False,
        bypass_url_encoding = True,
    )
    
    if (
        not os.path.exists(c_nocall.cache_file_name) or
        os.path.getsize(c_nocall.cache_file_name) == 0 or
        force_download
    ):
        
        c_login_1 = curl.Curl(
            urls.urls['msigdb']['login1'],
            cache = False,
            write_cache = False,
            follow = False,
            large = False,
            silent = True,
        )
        
        jsessionid = ''
        
        if hasattr(c_login_1, 'resp_headers'):
            
            for hdr in c_login_1.resp_headers:
                
                if hdr.lower().startswith(b'set-cookie'):
                    
                    jsessionid = hdr.split(b':')[1].split(b';')[0].strip()
                    jsessionid = jsessionid.decode('ascii')
                    _log('msigdb cookie obtained: `%s`.' % jsessionid)
                    break
        
        if not jsessionid:
            
            _log('msigdb: could not get cookie, returning empty list.')
            
            return {}
        
        req_headers = ['Cookie: %s' % jsessionid]
        
        c_login_2 = curl.Curl(
            urls.urls['msigdb']['login2'],
            cache = False,
            write_cache = False,
            large = False,
            silent = True,
            req_headers = req_headers,
            post = {
                'j_username': registered_email,
                'j_password': 'password',
            },
            follow = False,
            empty_attempt_again = False,
        )
        
        jsessionid_1 = ''
        
        if hasattr(c_login_2, 'resp_headers'):
            
            for hdr in c_login_2.resp_headers:
                
                if hdr.lower().startswith(b'set-cookie'):
                    
                    jsessionid_1 = hdr.split(b':')[1].split(b';')[0].strip()
                    jsessionid_1 = jsessionid_1.decode('ascii')
            
            _log(
                'msigdb: logged in with email `%s`, '
                'new cookie obtained: `%s`.' % (
                    registered_email,
                    jsessionid_1
                )
            )
        
        if not jsessionid_1:
            
            _log(
                'msigdb: could not log in with email `%s`, '
                'returning empty dict.' % registered_email
            )
            
            return {}
        
        req_headers_1 = ['Cookie: %s' % jsessionid_1]
    
    c = curl.Curl(
        url,
        req_headers = req_headers_1,
        silent = False,
        large = True,
        bypass_url_encoding = True,
        cache = not force_download,
    )
    
    result = {}
    
    for gset in c.result:
        
        gset = gset.strip().split('\t')
        
        result[gset[0]] = set(gset[2:])
    
    return result


def msigdb_download_collections(
        registered_email = None,
        only_collections = None,
        exclude = ('c5',),
        id_type = 'symbols',
    ):
    """
    Downloads all or some MSigDB gene set collections.
    Returns a dict of dicts where upper level keys are collections while
    lower level keys are geneset names and values are molecular identifiers.
    
    :arg str,NoneType registered_email:
        An email address registered at MSigDB. If `None` the `msigdb_email`
        from ``pypath.settings`` will be used.
    :arg set,NoneType only_collections:
        Limit the annotations only to these collections. For available
        collections e.g. ``{'h.all', 'c2cgp'}`` refer to the MSigDB webpage:
        http://software.broadinstitute.org/gsea/downloads.jsp#msigdb
    :arg tuple exclude:
        Exclude the collections having their name starting with any of the
        strings in this tuple. By default `c5` (Gene Ontology) is excluded.
    """
    
    MsigdbAnnotation = collections.namedtuple(
        'MsigdbAnnotation',
        [
            'collection',
            'geneset',
        ],
    )
    
    all_collections = {
        'hallmark': 'h.all',
        'positional': 'c1.all',
        'chemical_and_genetic_perturbations': 'c2.cgp',
        'biocarta_pathways': 'c2.cp.biocarta',
        'kegg_pathways': 'c2.cp.kegg',
        'pid_pathways': 'c2.cp.pid',
        'reactome_pathways': 'c2.cp.reactome',
        'mirna_targets': 'c3.mir',
        'tf_targets': 'c3.tft',
        'cancer_gene_neighborhoods': 'c4.cgn',
        'cancer_modules': 'c4.cm',
        'go_biological_process': 'c5.bp',
        'go_molecular_function': 'c5.mf',
        'go_cellular_component': 'c5.cc',
        'oncogenic_signatures': 'c6.all',
        'immunologic_signatures': 'c7.all',
    }
    
    collection_data = {}
    
    for collection, label in iteritems(all_collections):
        
        if (
            (
                only_collections and
                label not in only_collections
            ) or
            any(label.startswith(ex) for ex in exclude)
        ):
            
            continue
        
        _log(
            'MSigDB: downloading collection `%s` (%s).' % (collection, label)
        )
        
        collection_data[(collection, label)] = (
            msigdb_download(
                registered_email = registered_email,
                collection = label,
                id_type = id_type,
            )
        )
    
    return collection_data


def msigdb_annotations(
        registered_email = None,
        only_collections = None,
        exclude = ('c5',),
    ):
    """
    Downloads all or some MSigDB gene set collections and processes them
    to an annotation type dictionary.
    
    :arg str,NoneType registered_email:
        An email address registered at MSigDB. If `None` the `msigdb_email`
        from ``pypath.settings`` will be used.
    :arg set,NoneType only_collections:
        Limit the annotations only to these collections. For available
        collections e.g. ``{'h.all', 'c2cgp'}`` refer to the MSigDB webpage:
        http://software.broadinstitute.org/gsea/downloads.jsp#msigdb
    :arg tuple exclude:
        Exclude the collections having their name starting with any of the
        strings in this tuple. By default `c5` (Gene Ontology) is excluded.
     
    """
    
    MsigdbAnnotation = collections.namedtuple(
        'MsigdbAnnotation',
        [
            'collection',
            'geneset',
        ],
    )
    
    
    annotations = collections.defaultdict(set)
    
    collection_data = msigdb_download_collections(
        registered_email = registered_email,
        only_collections = only_collections,
        exclude = exclude,
    )
    
    for (collection, label), genesets in iteritems(collection_data):
        
        for geneset, genesymbols in iteritems(genesets):
            
            this_annot = MsigdbAnnotation(
                collection = collection,
                geneset = geneset,
            )
            
            for uniprot in mapping.map_names(
                genesymbols,
                'genesymbol',
                'uniprot',
            ):
                
                annotations[uniprot].add(this_annot)
    
    return dict(annotations)


def lrdb_interactions():
    
    
    resource_names = {
        'reactome': 'Reactome',
        'fantom5': 'Fantom5',
        'IUPHAR': 'Guide2Pharma',
        'uniprot': 'UniProt',
    }
    
    def remove(lst, to_remove):
        
        to_remove = common.to_set(to_remove)
        
        return [
            it
            for it in lst
            if it not in to_remove
        ]
    
    
    LrdbRecord = collections.namedtuple(
        'LrdbRecord',
        [
            'ligand_genesymbol',
            'receptor_genesymbol',
            'sources',
            'references',
            'ligand_cells',
            'receptor_cells',
        ]
    )
    
    
    url = urls.urls['lrdb']['url']
    
    c = curl.Curl(url, silent = False, large = True)
    
    reader = csv.DictReader(c.result, delimiter = '\t')
    
    result = []
    
    for rec in reader:
        
        result.append(
            LrdbRecord(
                ligand_genesymbol = rec['ligand'],
                receptor_genesymbol = rec['receptor'],
                sources = [
                    resource_names[src] if src in resource_names else src
                    for src in
                    remove(
                        rec['source'].split(','),
                        {'literature', ''},
                    )
                ],
                references = remove(rec['PMIDs'].split(','), ''),
                ligand_cells = remove(rec['cells.L'].split(','), ''),
                receptor_cells = remove(rec['cells.R'].split(','), ''),
            )
        )
    
    return result


def lrdb_annotations():
    
    LrdbAnnotation = collections.namedtuple(
        'LrdbAnnotation',
        [
            'role',
            'cell_type',
        ],
    )
    
    
    result = collections.defaultdict(set)
    
    
    lrdb = lrdb_interactions()
    
    for rec in lrdb:
        
        for role in ('ligand', 'receptor'):
            
            uniprots = mapping.map_name(
                getattr(rec, '%s_genesymbol' % role),
                'genesymbol',
                'uniprot',
            )
            
            for uniprot in uniprots:
                
                cell_types = getattr(rec, '%s_cells' % role) or (None,)
                
                for cell_type in cell_types:
                    
                    cell_type = (
                        'T lymphocyte'
                            if cell_type == 'tymphocyte' else
                        cell_type.replace('cells', 'cell')
                            if cell_type else
                        None
                    )
                    
                    result[uniprot].add(
                        LrdbAnnotation(
                            role = role,
                            cell_type = cell_type,
                        )
                    )
    
    return dict(result)


def baccin2019_interactions(ncbi_tax_id = 9606):
    
    def id_translate(mouse_gs):
        
        uniprots = mapping.map_name(
            mouse_gs,
            'genesymbol',
            'uniprot',
            10090,
        )
        
        if ncbi_tax_id != 10090:
            
            uniprots = set(
                itertools.chain(*(
                    homology.translate(uniprot)
                    for uniprot in uniprots
                ))
            )
        
        return uniprots
    
    
    def raw_to_uniprots(raw):
        
        components = raw.split('&')
        
        return set(
            itertools.product(
                *(id_translate(comp) for comp in components)
            )
        )
    
    
    def get_partners(components, sources, references):
        
        return {
            (
                comp[0]
                    if len(comp) == 1 else
                intera.Complex(
                    components = comp,
                    sources = sources,
                    references = references,
                )
            )
            for comp in components
        }
    
    
    Baccin2019Interaction = collections.namedtuple(
        'Baccin2019Interaction',
        [
            'ligand',
            'receptor',
            'correct',
            'ligand_location',
            'ligand_category',
            'resources',
            'references',
        ]
    )
    
    
    source_names = {
        'Baccin': 'Baccin2019',
        'Ramilowski': 'Ramilowski2015',
    }
    
    url = urls.urls['baccin2019']['url']
    
    c = curl.Curl(url, silent = False, large = True)
    
    data = inputs_common.read_xls(c.fileobj.name, sheet = 'SuppTable3')
    
    result = []
    
    if ncbi_tax_id != 10090:
        
        homology = homology_mod.ProteinHomology(
            target = ncbi_tax_id,
            source = 10090,
        )
    
    for rec in data[3:]:
        
        ligand_components = raw_to_uniprots(rec[1])
        
        if not ligand_components:
            
            continue
        
        receptor_components = raw_to_uniprots(rec[2])
        
        if not receptor_components:
            
            continue
        
        sources = {'Baccin2019', rec[3].strip()}
        sources = {
            source_names[s] if s in source_names else s
            for s in sources
        }
        
        references = {
            _ref for _ref in
            (
                ref.strip().replace('.0', '')
                for ref in rec[7].split(',')
            )
            if _ref.isdigit()
        }
        
        ligands = get_partners(ligand_components, sources, references)
        receptors = get_partners(receptor_components, sources, references)
        
        for ligand, receptor in itertools.product(ligands, receptors):
            
            result.append(
                Baccin2019Interaction(
                    ligand = ligand,
                    receptor = receptor,
                    correct = rec[4].strip(),
                    ligand_location = rec[5].strip(),
                    ligand_category = rec[6].strip(),
                    resources = sources,
                    references = references,
                )
            )
    
    return result


def baccin2019_annotations(ncbi_tax_id = 9606):
    
    Baccin2019Annotation = collections.namedtuple(
        'Baccin2019Annotation',
        [
            'mainclass',
            'subclass',
        ]
    )
    
    
    ia = baccin2019_interactions(ncbi_tax_id = ncbi_tax_id)
    
    result = collections.defaultdict(set)
    
    for ia in ia:
        
        result[ia.ligand].add(
            Baccin2019Annotation(
                mainclass = 'ligand',
                subclass = ia.ligand_category,
            )
        )
        
        result[ia.receptor].add(
            Baccin2019Annotation(
                mainclass = 'receptor',
                subclass = (
                    '%sReceptor' % ia.ligand
                        if ia.ligand_category != 'Other' else
                    None
                ),
            )
        )
    
    return dict(result)
