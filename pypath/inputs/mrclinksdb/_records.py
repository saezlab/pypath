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
    'MrclinksdbReceptorRaw',
    'MrclinksdbInteraction',
    'MrclinksdbReceptorInteraction',
    'MrclinksdbMetaboliteCell',
    'MrclinksdbTransporterRaw',
    'MrclinksdbTransporterInteraction',
]


MrclinksdbReceptorRaw = collections.namedtuple(
    'MrclinksdbReceptorRaw',
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

MrclinksdbReceptorInteraction = collections.namedtuple(
    'MrclinksdbReceptorInteraction',
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

# Backward-compatible aliases for the renamed receptor records.
MrclinksdbRaw = MrclinksdbReceptorRaw
MrclinksdbInteraction = MrclinksdbReceptorInteraction


MrclinksdbTransporterRaw = collections.namedtuple(
    'MrclinksdbTransporterRaw',
    [
        'hmdb_id',
        'metabolite_name',
        'hmdbp_id',        # HMDB protein ID; empty string for mouse (column absent)
        'enzyme_name',
        'gene_id',
        'gene_name',
        'uniprot_id',
        'type_',
    ],
)

MrclinksdbTransporterInteraction = collections.namedtuple(
    'MrclinksdbTransporterInteraction',
    [
        'hmdb',
        'metabolite_name',
        'hmdbp_id',
        'enzyme_name',
        'transporter_entrez',
        'transporter_genesymbol',
        'transporter_uniprot',
        'transporter_location',
    ],
)
