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

import collections

__all__ = [
    'SwisslipidsEvidence',
    'SwisslipidsLipid',
    'SwisslipidsTissue',
]


SwisslipidsLipid = collections.namedtuple(
    'SwisslipidsLipid',
    [
        'id',
        'level',
        'name',
        'abbreviation',
        'synonyms',
        'lipid_class',
        'parent',
        'components',
        'smiles',
        'inchi',
        'inchikey',
        'formula',
        'charge',
        'exact_mass',
        'chebi',
        'lipidmaps',
        'hmdb',
        'metanetx',
        'pmids',
    ]
)


SwisslipidsTissue = collections.namedtuple(
    'SwisslipidsTissue',
    [
        'lipid_id',
        'lipid_name',
        'ncbi_tax_id',
        'taxon_name',
        'tissue_uberon',
        'tissue_name',
        'evidence',
    ]
)


SwisslipidsEvidence = collections.namedtuple(
    'SwisslipidsEvidence',
    [
        'eco',
        'eco_name',
        'pmid',
    ]
)
