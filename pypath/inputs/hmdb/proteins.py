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

if TYPE_CHECKING:

    import pandas as pd

import pypath.inputs.hmdb.common as hmdb_common


def iter():
    """
    Itertate over protein records from HMDB.
    """

    return hmdb_common.iter(dataset = 'proteins')


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
            case see `pypath.inputs.hmdb.PROTEINS_SCHEMA` to see the full
            schema, and create your own that is restricted to your fields of
            interest.
        head:
            Process the first N records only. Useful for peeking into
            the data.

    Returns:
        A list of dicts, each dict containing the data extracted from
        one protein record. These dicts might be deeply nested, and
        the structure depends on the schema used.
    """

    return hmdb_common.raw(
        dataset = 'proteins',
        schema = schema,
        head = head,
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
           df = hmdb.proteins_table(
               'accession',
               'smiles',
               'state',
               head = 10,
           )
           df = hmdb.proteins_table(
               'accession',
               ('synonyms', '@'),
               head = 10,
           )
           df = hmdb.proteins_table(
               'accession',
               ('taxonomy', ('class', 'substituents')),
               head = 10,
           )

    """

    return hmdb_common.table(
        *fields,
        dataset = 'proteins',
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
        dataset = 'proteins',
        return_df = return_df,
        head = head,
    )
