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
import pypath.utils.taxonomy as taxonomy


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
        urls.urls['corum']['url_rescued'],
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
