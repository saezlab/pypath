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

import itertools

import pandas as pd

import pypath.inputs.hmdb.common as hmdb_common
from pypath.inputs.hmdb.schema.common import XMLNS, Field
from pypath.inputs.hmdb.schema.metabolites import SCHEMA
from pypath.inputs.hmdb import _log
import pypath.share.common as common
import pypath.formats.xml as xml


def hmdb_iter_metabolites():
    """
    Itertate over metabolite records from HMDB.
    """

    for _, element in hmdb_common.hmdb_xml('metabolites'):

        yield element


def hmdb_metabolites_raw(schema: dict = None, head: int | None = None) -> list:
    """
    Parse metabolite data from HMDB.

    Args:
        schema:
            The schema defines the fields to be parsed. By default, a schema
            covering nearly all fields in the HMDB metabolites XML is used.
            Likely you don't need all the fields within one task: in this
            case see `pypath.inputs.hmdb.METABOLITE_SCHEMA` to see the full
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

    result = []
    schema = schema or SCHEMA

    for i, met in enumerate(hmdb_iter_metabolites()):

        yield xml.fetch(met, schema, XMLNS)

        met.clear(keep_tail = True)

        if i == head:

            break

    return result


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

    for key in (id_type, f'{id_type}_id'):

        if key in SCHEMA:

            return key
