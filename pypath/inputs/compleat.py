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

import csv

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.internals.intera as intera
import pypath.utils.mapping as mapping


def compleat_raw():
    """
    Raw protein complex data from the Compleat database.
    """

    url = urls.urls['compleat']['rescued']
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

    return tab


def compleat_complexes(predicted = True):
    """
    Retrieves and processes protein complexes from the Compleat database.
    """

    raw = compleat_raw()
    complexes = {}

    for rec in raw:

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
