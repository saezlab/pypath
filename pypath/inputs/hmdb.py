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

from typing import Any, Literal

import itertools

from lxml import etree
import pandas as pd

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.common as common
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


def hmdb_metabolites_raw(schema: dict = None, head: int | None = None) -> list:
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

        yield xml.fetch(met, schema, XMLNS)

        met.clear(keep_tail = True)

        if i == head:

            break

    return result


class Field:

    def __init__(self, name, *definition):

        if isinstance(name, Field):

            name, definition = name.name, name.definition

        self.name = name
        self.d = definition or (name,)


    def process(self, record) -> tuple[tuple[Any]]:

        value = ((record,),)
        name = common.to_tuple(self.name)

        for d in self.d:

            if d == '@':

                value = tuple(
                    v[0]
                        if isinstance(v[0], (list, tuple)) else v
                    for v in value
                )

            elif d == '*' or isinstance(d, tuple):

                if d == '*':

                    d = tuple(sorted(value[0].keys()))

                if len(value) > 1:

                    _log('List of dicts content encountered.')

                value = value[0]
                value = tuple(tuple(v.get(k) for v in value) for k in d)
                name = tuple(f'{n}__{k}' for n in name for k in d)

            elif isinstance(d, str):

                value = tuple(tuple(v.get(d) for v in vv) for vv in value)

        return name, value


    def __str__(self):

        return self.name


def hmdb_table(
        *fields: str | tuple,
        head: int | None = None,
        **named_fields: str | tuple,
    ) -> pd.DataFrame:
    """
    Parse various simple and nested array fields from HMDB into data frame.

    Args:
        fields:
            Fields to include in the data frame. These must be keys in the
            schema, and will be also used as column names. Alternatively,
            tuples of sequetial processing steps can be provided: strings
            will be used as keys in nested dicts, tuples will be used as
            multiple keys in dicts, each yielding a separate column, the
            special symbol "*" means all keys in the sub-dict, while "@"
            means expand arrays into multiple rows. Be careful with this
            latter option because it is applied in a combinatorial way, i.e.
            in case of expanding an array to 5 rowns, and another one to 7
            rows results already 35 rows from a single record. This might
            result excessive memory use and processing time.
        named_fields:
            Same as `fields`, but the column name can be different from the
            top level key: argument names will be used as column names,
            values will be used as processing steps.
        head:
            Process the first N records only. Useful for peeking into
            the data.

    Examples:

        ..code-block:: python

           from pypath.inputs import hmdb
           df = hmdb.hmdb_table('accession', 'smiles', 'state', head = 10)
           df = hmdb.hmdb_table('accession', ('synonyms', '@'), head = 10)
           df = hmdb.hmdb_table(
               'accession',
               ('taxonomy', ('class', 'substituents')),
               head = 10,
           )

    """

    fields = [Field(d[0], *d) for d in (common.to_tuple(f) for f in fields)]
    fields.extend(Field(n, *f) for n, f in named_fields.items())
    keys = [f.d[0] for f in fields]
    schema = {k: SCHEMA[k] for k in keys}

    columns = []
    data = []
    result = pd.DataFrame()

    for i, met in enumerate(hmdb_metabolites_raw(schema, head = head)):

        cols, rows = zip(*(f.process(met) for f in fields))
        columns = columns or list(itertools.chain(*cols))
        data.extend(itertools.product(*itertools.chain(*rows)))

        if (i + 1) % 100 == 0:

            df = pd.DataFrame.from_records(data, columns = columns)
            result = pd.concat((result, df))
            data = []

    df = pd.DataFrame.from_records(data, columns = columns)
    result = pd.concat((result, df))

    return result


def hmdb_mapping(
        id_type_a: str,
        id_type_b: str,
        return_df: bool = False,
        head: int | None = None,
    ) -> dict[str, set[str]] | pd.DataFrame:
    """
    ID translation input from HMDB.

    Note: you can use this function for purposes other than ID translation
    tables, e.g. you can collect the molecular weights, pathways, etc. Though
    not all these options are guaranteed to work or result a meaningful
    output.

    Args:
        id_type_a:
            An identifier type, see the `ID_FIELDS` set in this module,
            and keys in the `SIMPLE_FIELDS` dict.
        id_type_b:
            Another identifier type, same options as for `id_type_a`.
        return_df:
            Return a data frame instead of dict of sets.
        head:
            Process the first N records only. Useful for peeking into
            the data.

    Return:
        Translation data between two types of identifiers.
    """

    fields = {
        'id_a': (_id_type(id_type_a), '@'),
        'id_b': (_id_type(id_type_b), '@'),
    }

    df = hmdb_table(**fields, head = head)

    if return_df:

        return df

    else:

        return df.groupby('id_a')['id_b'].apply(set).to_dict()


def _id_type(id_type: str) -> str:
    """
    Field name from ID type.
    """

    return SCHEMA.get('id_type', SCHEMA.get(f'{id_type}_id'))
