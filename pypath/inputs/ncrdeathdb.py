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
import collections

import pypath.resources.urls as urls
import pypath.share.curl as curl


def ncrdeathdb_interactions():

    NcrdeathdbInteraction = collections.namedtuple(
        'NcrdeathdbInteraction',
        (
            'ncrna',
            'target_gene',
            'ncrna_type',
            'pathway',
            'effect',
            'pmid',
            'organism',
        ),
    )

    url = urls.urls['ncrdeathdb']['url_rescued']
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

            result.append(
                NcrdeathdbInteraction(
                    ncrna = rna_id,
                    target_gene = protein_id,
                    ncrna_type = typ,
                    pathway = rec['Pathway'].strip(),
                    effect = rec['Action_Mode'].strip() or None,
                    pmid = rec['PMID'].strip(),
                    organism = int(rec['tax_id'].strip()),
                )
            )

    return result
