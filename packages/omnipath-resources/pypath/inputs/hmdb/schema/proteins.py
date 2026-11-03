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

import lxml.etree as etree

import pypath.formats.xml as xml
from pypath.inputs.hmdb.schema.common import (
    XMLNS,
    SIMPLE_FIELDS,
    ARRAY_FIELDS,
    PATHWAYS,
    _ref_chunk,
)


SIMPLE_FIELDS = SIMPLE_FIELDS.copy()
SIMPLE_FIELDS.update({
    'general_function',
    'specific_function',
    'protein_type',
})

ID_FIELDS = {
    'gene_name',
    'genebank_protein_id',
    'uniprot_id',
    'uniprot_name',
    'hgnc_id',
    'genecard_id',
    'geneatlas_id',
}

SIMPLE_FIELDS.update(ID_FIELDS)
SIMPLE_FIELDS = xml._simple_fields(SIMPLE_FIELDS)

ARRAY_FIELDS = ARRAY_FIELDS.copy()
ARRAY_FIELDS.update({
    ('subcellular_locations', 'subcellular_location'),
    ('pdb_ids', 'pdb_id'),
})
ARRAY_FIELDS = xml._array_fields(ARRAY_FIELDS)

GENE_PROPERTIES = (
    'gene_properties',
    xml._simple_fields((
        'chromosome_location',
        'locus',
        'gene_sequence',
    ))
)

PROTEIN_PROPERTIES = (
    'protein_properties',
    dict(
        **xml._simple_fields((
            'residue_number',
            'molecular_weight',
            'theoretical_pi',
            'polypeptide_sequence',
        )),
        **xml._array_fields((
            ('transmembrane_regions', 'region'),
            ('signal_regions', 'region'),
        )),
    ),
)

PFAMS = (
    'pfams',
    ('pfam', 'findall'),
    {
        'name',
        'pfam_id',
    },
)

METABOLITE_ASSOCIATIONS = (
    'metabolite_associations',
    ('metabolite', 'findall'),
    {
        'accession',
        'name',
    },
)

GO_CLASSIFICATIONS = (
    'go_classifications',
    ('go_class', 'findall'),
    {'category', 'description', 'go_id'},
)

GENERAL_REFERENCES = _ref_chunk('general_references')

METABOLITE_REFERENCES = (
    'metabolite_references',
    ('metabolite_reference', 'findall'),
    {
        'metabolite': (
            'metabolite',
            {'name', 'accession'},
        ),
        'reference': (
            'reference',
            'pubmed_id',
            None,
        ),
    },
)

SCHEMA = {
    'gene_properties': GENE_PROPERTIES,
    'protein_properties': PROTEIN_PROPERTIES,
    'pfams': PFAMS,
    'metabolite_associations': METABOLITE_ASSOCIATIONS,
    'go_classifications': GO_CLASSIFICATIONS,
    'pathways': PATHWAYS,
    'general_references': GENERAL_REFERENCES,
    'metabolite_references': METABOLITE_REFERENCES,
}
SCHEMA.update(SIMPLE_FIELDS)
SCHEMA.update(ARRAY_FIELDS)
