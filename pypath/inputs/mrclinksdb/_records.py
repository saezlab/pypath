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

import collections

__all__ = [
    'MrclinksdbRaw',
    'MrclinksdbInteraction',
    'MrclinksdbMetaboliteCell',
]


MrclinksdbRaw = collections.namedtuple(
    'MrclinksdbRaw',
    [
        'mrid',
        'hmdb_id',
        'metabolite_name',
        'pubchem_cid_sid',
        'molecular_formula',
        'kingdom',
        'super_class_',
        'class_',
        'canonical_smiles',
        'receptor_gene_id',
        'receptor_uniprot_id',
        'receptor_symbol',
        'protein_name',
        'pmid',
        'other_db_',
    ]
)

MrclinksdbInteraction = collections.namedtuple(
    'MrclinksdbInteraction',
    [
        'mrid',
        'hmdb',
        'name',
        'pubchem',
        'pubchem_sid',
        'formula',
        'compound_kingdom',
        'compound_superclass',
        'compound_class',
        'smiles',
        'receptor_entrez',
        'receptor_uniprot',
        'receptor_genesymbol',
        'pmids',
        'resource',
        'receptor_location',
    ],
)

MrclinksdbMetaboliteCell = collections.namedtuple(
    'MrclinksdbMetaboliteCell',
    [
        'hmdb',
        'metabolite',
        'interaction',
        'cell_type',
        'experimental_subject',
        'disease',
        'effect',
        'pmids',
    ],
)
