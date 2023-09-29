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

import bs4

import pypath.share.curl as curl
import pypath.share.progress as progress
import pypath.resources.urls as urls
import pypath.utils.mapping as mapping
import pypath.inputs.uniprot_db as uniprot_db


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

    soup = bs4.BeautifulSoup(c.result, features = 'xml')
    uniprots = uniprot_db.get_db(organism = organism, swissprot = None)
    prg = progress.Progress(
        len(soup.find_all('record')),
        'Processing data from LMPID',
        21
    )

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
    ``pypath.inputs.lmpid.load_lmpid`` to list of interactions.
    """

    data = load_lmpid(organism = organism)

    return [[l['prey'], l['bait'], ';'.join(l['refs'])] for l in data]


def lmpid_dmi(organism = 9606):
    """
    Converts list of domain-motif interactions supplied by
    ``pypath.inputs.lmpid.load_lmpid`` to list of
    ``pypath.intera.DomainMotif`` objects.
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
