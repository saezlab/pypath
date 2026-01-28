#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright 2014-2024
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

"""
Named tuple definitions for Metabolic Atlas data structures.
"""

from __future__ import annotations

import collections

__all__ = [
    'MetatlasModel',
    'MetatlasGem',
    'MetatlasReaction',
    'MetatlasMetabolite',
    'MetatlasGene',
]


MetatlasModel = collections.namedtuple(
    'MetatlasModel',
    [
        'id',
        'name',
        'short_name',
        'organism',
        'tissue',
        'cell_type',
        'condition',
        'year',
        'reaction_count',
        'metabolite_count',
        'gene_count',
        'files',
    ],
)


MetatlasGem = collections.namedtuple(
    'MetatlasGem',
    [
        'model_name',
        'git_host',
        'git_repo',
        'git_url',
        'latest_version',
        'commits',
        'contributors',
    ],
)


MetatlasReaction = collections.namedtuple(
    'MetatlasReaction',
    [
        'id',
        'kegg',
        'bigg',
        'metanetx',
        'rhea',
        'rhea_master',
        'reactome',
        'spontaneous',
    ],
)


MetatlasMetabolite = collections.namedtuple(
    'MetatlasMetabolite',
    [
        'id',
        'id_no_compartment',
        'bigg',
        'kegg',
        'hmdb',
        'chebi',
        'pubchem',
        'lipidmaps',
        'metanetx',
    ],
)


MetatlasGene = collections.namedtuple(
    'MetatlasGene',
    [
        'id',
        'ensembl_transcript',
        'ensembl_protein',
        'uniprot',
        'symbol',
        'entrez',
        'name',
        'aliases',
        'compartments',
    ],
)
