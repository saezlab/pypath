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


ONTOLOGY_FIELDS = xml._simple_fields((
    'term',
    'definition',
    'parent_id',
    'level',
    'type',
))


def _ontology(descendant: etree._Element) -> dict:

    path = (
        dict(
            **ONTOLOGY_FIELDS,
            **{
                'descendants': (
                    'descendants',
                    ('descendant', 'findall'),
                    _ontology,
                ),
            },
        ),
    )

    return xml.fetch(descendant, path, namespaces = XMLNS)


ONTOLOGY = (
    'ontology',
    ('root', 'findall'),
    _ontology,
)

SIMPLE_FIELDS = SIMPLE_FIELDS.copy()
SIMPLE_FIELDS.update({
    'accession',
    'name',
    'version',
    'status',
    'description',
    'chemical_formula',
    'average_molecular_weight',
    'monisotopic_molecular_weight',
    'iupac_name',
    'traditional_iupac',
    'cas_registry_number',
    'smiles',
    'inchi',
    'inchikey',
    'state',
    'synthesis_reference',
})

ID_FIELDS = {
    'chemspider',
    'drugbank',
    'foodb',
    'pubchem_compound',
    'pdb',
    'chebi',
    'phenol_explorer_compound',
    'knapsack',
    'kegg',
    'biocyc',
    'bigg',
    'wikipedia',
    'metlin',
    'vmh',
    'fbonto',
}

SIMPLE_FIELDS.update(xml._simple_fields(f'{f}_id' for f in ID_FIELDS))
ARRAY_FIELDS = ARRAY_FIELDS.copy()

TAXONOMY = (
    'taxonomy',
    dict(
        **xml._simple_fields((
            'description',
            'direct_parent',
            'kingdom',
            'class',
            'sub_class',
            'molecular_framework',
        )),
        **xml._array_fields((
            ('alternative_parents', 'alternative_parent'),
            ('substituents', 'substituent'),
        )),
    ),
)

EXPERIMENTAL_PROPERTIES = (
    'experimental_properties',
    ('property', 'findall'),
    {'kind', 'value', 'source'},
)

PREDICTED_PROPERTIES = (
    'predicted_properties',
    ('property', 'findall'),
    {'kind', 'value', 'source'},
)

SPECTRA = (
    'spectra',
    ('spectrum', 'findall'),
    {'type', 'spectrum_id'},
)

BIOLOGICAL_PROPERTIES = (
    'biological_properties',
    {
        'cellular_locations': (
            'cellular_locations',
            ('cellular', 'findall'),
            None,
        ),
        'biospecimen_locations': (
            'biospecimen_locations',
            ('biospecimen', 'findall'),
            None,
        ),
        'tissue_locations': (
            'tissue_locations',
            ('tissue', 'findall'),
            None,
        ),
        'pathways': PATHWAYS,
    },
)

NORMAL_CONCENTRATIONS = (
    'normal_concentrations',
    ('concentration', 'findall'),
    dict(
        **xml._simple_fields((
            'biospecimen',
            'concentration_value',
            'concentration_units',
            'subject_age',
            'subject_sex',
            'subject_condition',
        )),
        **{'references': _ref_chunk()},
    ),
)

ABNORMAL_CONCENTRATIONS = (
    'abnormal_concentrations',
    ('concentration', 'findall'),
    dict(
        **xml._simple_fields((
            'biospecimen',
            'concentration_value',
            'concentration_units',
            'patient_age',
            'patient_sex',
            'patient_information',
        )),
        **{'references': _ref_chunk()},
    ),
)

DISEASES = (
    'diseases',
    ('disease', 'findall'),
    dict(
        **xml._simple_fields((
            'name',
            'omim_id',
        )),
        **{'references': _ref_chunk()},
    ),
)

GENERAL_REFERENCES = _ref_chunk('general_references')

PROTEIN_ASSOCIATIONS = (
    'protein_associations',
    ('protein', 'findall'),
    {
        'protein_accession',
        'name',
        'uniprot_id',
        'gene_name',
        'protein_type',
    },
)

SCHEMA = {
    'taxonomy': TAXONOMY,
    'spectra': SPECTRA,
    'biological_properties': BIOLOGICAL_PROPERTIES,
    'experimental_properties': EXPERIMENTAL_PROPERTIES,
    'predicted_properties': PREDICTED_PROPERTIES,
    'normal_concentrations': NORMAL_CONCENTRATIONS,
    'abnormal_concentrations': ABNORMAL_CONCENTRATIONS,
    'diseases': DISEASES,
    'general_references': GENERAL_REFERENCES,
    'protein_associations': PROTEIN_ASSOCIATIONS,
    'ontology': ONTOLOGY,
}
SCHEMA.update(xml._simple_fields(SIMPLE_FIELDS))
SCHEMA.update(xml._array_fields(ARRAY_FIELDS))
