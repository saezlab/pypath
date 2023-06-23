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
import csv
import collections

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.share.common as common
import pypath.share.settings as settings


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
    c = curl.Curl(
        url,
        large = True,
        silent = False,
        files_needed = files,
        req_headers = [settings.get('user_agent')],
        alpn = False,
    )

    evidences = collections.defaultdict(list)

    for rec in csv.DictReader(c.files_multipart['evidences.csv']):

        evidences[rec['ID']].append(rec)

    records = list(csv.DictReader(c.files_multipart['export.csv']))

    return records, evidences


def protmapper_enzyme_substrate(
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
        'signor': 'SIGNOR',
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
            databases[source] if source in databases else source
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

    return protmapper_enzyme_substrate(interactions = True, **kwargs)
