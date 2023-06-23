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

from past.builtins import xrange, range

import re
import collections

import bs4

import pypath.share.curl as curl
import pypath.resources.urls as urls


def elm_domains():

    url = urls.urls['ielm_domains']['url']
    c = curl.Curl(url, silent = False)
    data = c.result
    soup = bs4.BeautifulSoup(data, 'html.parser')
    tbl = soup.find('table').find_all('td')
    rows = [tbl[x:x + 4] for x in xrange(0, len(tbl), 4)]
    result = {}

    for r in rows:

        uniprot = r[1].text
        motif = r[0].text

        if uniprot not in result:

            result[uniprot] = {}

        if motif not in result[uniprot]:

            result[uniprot][motif] = []

        result[uniprot][motif].append((r[2].text, r[3].text))

    return result


def elm_classes():

    ELMClass = collections.namedtuple(
        'ELMClass',
        (
            'accession',
            'identifier',
            'functional_name',
            'description',
            'regex',
            'probability',
            'n_instances',
            'n_pdb',
        ),
    )

    url = urls.urls['elm_class']['url']
    c = curl.Curl(url, silent = False)
    data = c.result
    data = [
        ELMClass(*x.split('\t'))
        for x in data.replace('"', '').split('\n')[6:]
        if len(x) > 0
    ]

    return dict(zip(
        (x[1] for x in data),
        data,
    ))


def elm_instances():

    ELMInstance = collections.namedtuple(
        'ELMInstance',
        (
            'accession',
            'type',
            'identifier',
            'uniprot_id',
            'uniprot',
            'synonyms',
            'start',
            'end',
            'references',
            'methods',
            'logic',
            'pdb',
            'organism',
        ),
    )

    url = urls.urls['elm_inst']['url']
    c = curl.Curl(url, silent = False, slow = True)
    data = c.result
    data = data.replace('"', '').split('\n')
    data = [
        x.split('\t')
        for x in data[6:]
    ]

    return [
        ELMInstance(*x)
        for x in data
        if len(x) == 13
    ]


def elm_interactions():
    """
    Downlods manually curated interactions from ELM.
    This is the gold standard set of ELM.
    """

    def number_or_none(value, typ = int):

        return typ(value) if value != 'None' else None

    # UniProt ID with isoform e.g. O14754-1
    reupi = re.compile(r'([\w]{6,10})(?:-([0-9]{1,2}))?')
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
    c = curl.Curl(url, silent = False, slow = True)
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
