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
Access the Human Metabolome Database.
"""

from typing import Literal

import itertools

import pandas as pd

import pypath.inputs.hmdb.xml as hmdb_xml
from pypath.inputs.hmdb.schema.common import XMLNS, Field
from pypath.inputs.hmdb.schema.metabolites import SCHEMA as METABOLITES_SCHEMA
from pypath.inputs.hmdb.schema.proteins import SCHEMA as PROTEINS_SCHEMA
import pypath.share.common as common
import pypath.formats.xml as xml


def iter(dataset: Literal['metabolites', 'proteins']):
    """
    Itertate over records from HMDB.
    """

    for _, element in hmdb_xml.hmdb_xml(dataset):

        yield element


def raw(
        dataset: Literal['metabolites', 'proteins'],
        schema: dict = None,
        head: int | None = None,
    ) -> list[dict]:
    """
    Parse data from HMDB XML files.

    Args:
        dataset:
            Either "metabolites" or "proteins".
        schema:
            The schema defines the fields to be parsed. By default, a schema
            covering nearly all fields in the HMDB XML is used.
            Likely you don't need all the fields within one task: in this
            case see `pypath.inputs.hmdb.METABOLITE_SCHEMA` or
            `pypath.inputs.hmdb.PROTEINS_SCHEMA` to see the full
            schema, and create your own that is restricted to your fields of
            interest.
        head:
            Process the first N records only. Useful for peeking into
            the data.

    Returns:
        A list of dicts, each dict containing the data extracted from
        one record. These dicts might be deeply nested, and the structure
        depends on the schema used.
    """

    result = []
    schema = schema or globals()[f'{dataset.upper()}_SCHEMA']

    for i, record in enumerate(iter(dataset)):

        yield xml.fetch(record, schema, XMLNS)

        record.clear(keep_tail = True)

        if i == head:

            break

    return result


def processed(
        *fields: str | tuple,
        dataset: Literal['metabolites', 'proteins'],
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
    """

    SCHEMA = globals()[f'{dataset.upper()}_SCHEMA']
    fields = [Field(d[0], *d) for d in (common.to_tuple(f) for f in fields)]
    fields.extend(Field(n, *f) for n, f in named_fields.items())
    keys = [f.d[0] for f in fields]
    schema = {k: SCHEMA[k] for k in keys}
    columns = []

    if not fields:

        raise ValueError('At least one field must be provided.')

    for record in raw(dataset, schema, head = head):

        cols, rows = zip(*(f.process(record) for f in fields))
        columns = columns or list(itertools.chain(*cols))

        for record in itertools.product(*itertools.chain(*rows)):

            yield record, columns


def table(
        *fields: str | tuple,
        dataset: Literal['metabolites', 'proteins'],
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
    """

    data = []
    result = pd.DataFrame()

    for i, (rec, cols) in enumerate(processed(
        *fields,
        dataset = dataset,
        head = head,
        **named_fields
    )):

        data.append(rec)

        if (i + 1) % 1000 == 0:

            df = pd.DataFrame.from_records(data, columns = cols)
            result = pd.concat((result, df))
            data = []

    df = pd.DataFrame.from_records(data, columns = cols)
    result = pd.concat((result, df))

    return result


def mapping(
        id_type_a: str,
        id_type_b: str,
        dataset: Literal['metabolites', 'proteins'],
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

    df = table(**fields, dataset = dataset, head = head)

    if return_df:

        return df

    else:

        return df.groupby('id_a')['id_b'].apply(set).to_dict()


def _id_type(id_type: str) -> str:
    """
    Field name from ID type.
    """

    for key in (id_type, f'{id_type}_id'):

        if key in METABOLITES_SCHEMA:

            return key
