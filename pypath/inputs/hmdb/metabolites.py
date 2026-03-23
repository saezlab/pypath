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

from __future__ import annotations

"""
Preprocessed metabolite data from the Human Metabolome Database.
"""

from typing import TYPE_CHECKING
from collections import namedtuple

if TYPE_CHECKING:

    import pandas as pd

import pypath.inputs.hmdb.common as hmdb_common




def iter():
    """
    Itertate over metabolite records from HMDB.
    """

    return hmdb_common.iter(dataset = 'metabolites')


def raw(
        schema: dict = None,
        head: int | None = None,
    ) -> list[dict]:
    """
    Parse metabolite data from HMDB.

    Args:
        schema:
            The schema defines the fields to be parsed. By default, a schema
            covering nearly all fields in the HMDB metabolites XML is used.
            Likely you don't need all the fields within one task: in this
            case see `pypath.inputs.hmdb.METABOLITES_SCHEMA` to see the full
            schema, and create your own that is restricted to your fields of
            interest.
        head:
            Process the first N records only. Useful for peeking into
            the data.

    Returns:
        A list of dicts, each dict containing the data extracted from
        one metabolite record. These dicts might be deeply nested, and
        the structure depends on the schema used.
    """

    return hmdb_common.raw(
        dataset = 'metabolites',
        schema = schema,
        head = head,
    )


def processed(
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
           mets = list(hmdb.metabolites_processed(
               'accession',
               'smiles',
               'state',
               head = 10,
           ))
           mets = list(hmdb.metabolites_processed(
               'accession',
               ('synonyms', '@'),
               head = 10,
           ))
           mets = list(hmdb.metabolites_processed(
               'accession',
               ('taxonomy', ('class', 'substituents')),
               head = 10,
           ))

    """

    return hmdb_common.processed(
        *fields,
        dataset = 'metabolites',
        head = head,
        **named_fields,
    )


def table(
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
           df = hmdb.metabolites_table(
               'accession',
               'smiles',
               'state',
               head = 10,
           )
           df = hmdb.metabolites_table(
               'accession',
               ('synonyms', '@'),
               head = 10,
           )
           df = hmdb.metabolites_table(
               'accession',
               ('taxonomy', ('class', 'substituents')),
               head = 10,
           )

    """

    return hmdb_common.table(
        *fields,
        dataset = 'metabolites',
        head = head,
        **named_fields,
    )


def mapping(
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

    return hmdb_common.mapping(
        id_type_a = id_type_a,
        id_type_b = id_type_b,
        dataset = 'metabolites',
        return_df = return_df,
        head = head,
    )


CompoundsForMetabo = namedtuple('CompoundsForMetabo', [
    'synonyms',
    'chebi_id',
    'accession',
    'pubchem_compound_id',
    'kegg_id',
    'drugbank_id',
    'cas_registry_number',
    'traditional_iupac',
    'iupac_name',
    'monisotopic_molecular_weight',
    'average_molecular_weight',
    'chemical_formula',
    'inchi',
    'inchikey',
    'smiles',
    'general_references',
])


def synonyms_chebi() -> dict[str, str]:
    """
    Compound name and synonym → ChEBI ID mapping from HMDB.

    Parses the HMDB metabolites XML once, collecting primary name,
    IUPAC name, traditional IUPAC name, and all synonyms for each
    compound that carries a ChEBI annotation.  Keys are lower-cased
    for case-insensitive lookup.

    Returns:
        Dict mapping lowercase name/synonym strings to ChEBI ID strings
        (e.g. ``'atp'`` → ``'CHEBI:30616'``).  When multiple compounds
        share a synonym, the first encountered mapping is kept.
    """

    from pypath.inputs.hmdb.schema.metabolites import SCHEMA

    _name_fields = ('name', 'iupac_name', 'traditional_iupac')
    sub_schema = {
        k: SCHEMA[k]
        for k in (*_name_fields, 'synonyms', 'chebi_id')
    }

    result: dict[str, str] = {}

    for rec in hmdb_common.raw('metabolites', schema = sub_schema):

        chebi = rec.get('chebi_id') or ''

        if not chebi:
            continue

        if not str(chebi).upper().startswith('CHEBI:'):
            chebi = f'CHEBI:{chebi}'

        for field in _name_fields:
            val = rec.get(field) or ''
            if val:
                result.setdefault(str(val).lower(), chebi)

        for syn in (rec.get('synonyms') or []):
            if syn:
                result.setdefault(str(syn).lower(), chebi)

    return result


def compounds_for_metabo():

    for rec in processed(
        'synonyms',
        'chebi_id',
        'accession',
        'pubchem_compound_id',
        'kegg_id',
        'drugbank_id',
        'cas_registry_number',
        'traditional_iupac',
        'iupac_name',
        'monisotopic_molecular_weight',
        'average_molecular_weight',
        'chemical_formula',
        'inchi',
        'inchikey',
        'smiles',
        'general_references',
    ):

        yield CompoundsForMetabo(*rec[0])
