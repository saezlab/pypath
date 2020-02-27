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

import re
import collections
import itertools

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.utils.mapping as mapping
import pypath.inputs.common as inputs_common


def rolland_hi_ii_14():
    """
    Loads the HI-II-14 unbiased interactome from the large scale screening
    of from Rolland 2014.
    Returns list of interactions.
    """
    url = urls.urls['hiii14']['url']
    c = curl.Curl(url, silent = False, large = True)
    xlsname = c.fileobj.name
    c.fileobj.close()
    tbl = inputs_common.read_xls(xlsname, sheet = '2G')

    for row in tbl[1:]:

        yield [c.split('.')[0] for c in row]


def vidal_hi_iii_old(fname):
    """
    Loads the HI-III  unbiased interactome from preliminary data of
    the next large scale screening of Vidal Lab.

    The data is accessible here:
        http://interactome.dfci.harvard.edu/H_sapiens/dload_trk.php
    You need to register and accept the license terms.

    Returns list of interactions.
    """

    f = curl.FileOpener(fname)

    return [l.strip().split('\t') for l in f.result][1:]


def hi_iii():
    """
    Loads the unbiased human interactome version III (HI-III).
    This is an unpublished data and its use is limited.
    Please check the conditions and licensing terms carefully at
    http://interactome.baderlab.org.
    """

    HiiiiInteraction = collections.namedtuple(
        'HiiiiInteraction',
        [
            'id_a',
            'id_b',
            'isoform_a',
            'isoform_b',
            'screens',
            'score',
        ]
    )


    rescore = re.compile(r'author score: ([\d\.]+)')
    rescreens = re.compile(r'Found in screens ([\d,]+)')

    url = urls.urls['hid']['hi-iii']
    post_data = {
        'form[request_dataset]': '2',
        'form[request_file_format]': 'psi',
    }
    c = curl.Curl(url, silent = False, large = True, post = post_data)

    for row in c.result:

        if not row.strip():

            continue

        id_a, id_b, rest = row.split(' ', maxsplit = 2)
        id_a, isoform_a = id_a.split('-') if '-' in id_a else (id_a, 1)
        id_b, isoform_b = id_b.split('-') if '-' in id_b else (id_b, 1)

        sc = rescore.search(rest)
        score = float(sc.groups()[0]) if sc else None
        screens = tuple(
            int(i) for i in rescreens.search(rest).groups()[0].split(',')
        )

        yield HiiiiInteraction(
            id_a = id_a[10:],
            id_b = id_b[10:],
            isoform_a = int(isoform_a),
            isoform_b = int(isoform_b),
            screens = screens,
            score = score,
        )


def lit_bm_13_interactions():
    """
    Downloads and processes Lit-BM-13 dataset, the 2013 version of the
    high confidence literature curated interactions from CCSB.
    Returns list of interactions.
    """

    LitBm13Interaction = collections.namedtuple(
        'LitBm13Interaction',
        [
            'entrez_a',
            'entrez_b',
            'genesymbol_a',
            'genesymbol_b',
        ]
    )

    url = urls.urls['hid']['lit-bm-13']
    c = curl.Curl(url, silent = False, large = True)

    _ = next(c.result)

    for row in c.result:

        row = row.strip().split('\t')

        yield LitBm13Interaction(
            entrez_a = row[0],
            entrez_b = row[2],
            genesymbol_a = row[1],
            genesymbol_b = row[3],
        )


def lit_bm_17_interactions():
    """
    Downloads and processes Lit-BM-13 dataset, the 2017 version of the
    high confidence literature curated interactions from CCSB.
    Returns list of interactions.
    """

    LitBm17Interaction = collections.namedtuple(
        'LitBm17Interaction',
        [
            'id_a',
            'id_b',
            'pubmed',
            'score',
        ]
    )

    url = urls.urls['hid']['lit-bm-17']
    c = curl.Curl(url, silent = False)
    data = c.result

    c = curl.Curl(url, silent = False, large = True)

    _ = next(c.result)

    for row in c.result:

        row = row.strip().split('\t')

        id_a = row[0][10:]
        id_b = row[1][10:]

        pubmed = row[8][7:]
        score = float(row[14][13:])

        yield LitBm17Interaction(
            id_a = id_a,
            id_b = id_b,
            pubmed = pubmed,
            score = score,
        )


def huri_interactions():

    return _huri_interactions(dataset = 'huri')


def yu2011_interactions():

    return _huri_interactions(dataset = 'yu-2011')


def hi_union_interactions():

    return _huri_interactions(dataset = 'hi-union')


def yang2016_interactions():

    return _huri_interactions(dataset = 'yang-2016')


def _huri_interactions(dataset):

    reuniprot = re.compile(r'[a-z]+:([\w\.]+)(?:-?([0-9]?))?')
    rescore = re.compile(r'author score: ([\.0-9]+)')

    HuriInteraction = collections.namedtuple(
        'HuriInteraction',
        [
            'uniprot_a',
            'uniprot_b',
            'isoform_a',
            'isoform_b',
            'score',
        ]
    )


    def _map_ids(_id):

        return mapping.map_name(
            _id,
            _id[:4].lower() if _id[:4] in {'ensp', 'enst'} else 'uniprot',
            'uniprot',
        )


    url = dataset if dataset.startswith('http') else urls.urls['hid'][dataset]
    c = curl.Curl(url, large = True, silent = False)

    for row in c.result:

        score = rescore.search(row)

        if score:

            score = float(score.groups()[0])

        row = row.split()

        if len(row) < 2:

            continue

        id_a, isoform_a = reuniprot.match(row[0]).groups()
        id_b, isoform_b = reuniprot.match(row[1]).groups()

        uniprots_a = _map_ids(id_a)
        uniprots_b = _map_ids(id_b)

        for uniprot_a, uniprot_b in itertools.product(uniprots_a, uniprots_b):

            yield HuriInteraction(
                uniprot_a = uniprot_a,
                uniprot_b = uniprot_b,
                isoform_a = int(isoform_a) if isoform_a else 1,
                isoform_b = int(isoform_b) if isoform_b else 1,
                score = score,
            )
