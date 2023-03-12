#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2023
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: Dénes Türei (turei.denes@gmail.com)
#           Nicolàs Palacio
#           Sebastian Lobentanzer
#           Erva Ulusoy
#           Olga Ivanova
#           Ahmet Rifaioglu
#           Sebastian Lobentanzer
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

from __future__ import annotations

"""
Access the Human Metabolome Database.
"""

from typing import Literal

from lxml import etree

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.session as session
import pypath.formats.xml as xml

_log = session.Logger(name = 'hmdb_input')._log


XMLNS = '{http://www.hmdb.ca}'


def hmdb_xml(dataset: Literal['metabolites']) -> etree.iterparse:
    """
    Download and open the XML file of a dataset from the HMDB.
    """

    url = urls.urls['hmdb'][dataset]
    c = curl.Curl(
        url,
        large = True,
        silent = False,
        slow = True,
        default_mode = 'rb',
    )

    return etree.iterparse(
        c.result['hmdb_metabolites.xml'],
        tag = f'{XMLNS}metabolite',
    )


def hmdb_iter_metabolites():
    """
    Itertate over metabolite records from HMDB.
    """

    for _, element in hmdb_xml('metabolites'):

        yield element


def _ref_chunk(container: str = 'references') -> tuple:

    return (
        container,
        ('reference', 'findall'),
        'pubmed_id',
        None,
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

SIMPLE_FIELDS = {
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
}

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

ARRAY_FIELDS = {
    ('secondary_accessions', 'accession'),
    ('synonyms', 'synonym'),
}

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
        'pathways': (
            'pathways',
            ('pathway', 'findall'),
            {'name', 'smpdb_id', 'kegg_map_id'},
        ),
    }
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


def hmdb_metabolites(schema: dict = None, head: int | None = None) -> list:
    """
    Parse metabolite data from HMDB.

    Args:
        schema:
            The schema defines the fields to be parsed. By default, a schema
            covering nearly all fields in the HMDB metabolites XML is used.
            Likely you don't need all the fields within one task: in this
            case see `pypath.inputs.hmdb.SCHEMA` to see the full schema,
            and create your own that is restricted to your fields of interst.
        head:
            Process the first N records only. Useful for peeking into
            the data.

    Returns:
        A list of dicts, each dict containing the data extracted from
        one metabolite record. These dicts might be deeply nested, and
        the structure depends on the schema used.
    """

    result = []
    schema = schema or SCHEMA

    for i, met in enumerate(hmdb_iter_metabolites()):

        result.append(xml.fetch(met, schema, XMLNS))

        if i == head:

            break

    return result
