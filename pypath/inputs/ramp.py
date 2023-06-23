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
Access the RaMP metabolomic pathway and metabolite database.
"""

from typing import IO, Literal, TYPE_CHECKING

if TYPE_CHECKING:

    import sqlite3

import os
import json
import pprint

import pandas as pd

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.session as session
import pypath.share.common as common
import pypath.formats.sqldump as sqldump

_log = session.Logger(name = 'ramp_input')._log


def _ramp_sqldump() -> IO:
    """
    Download the RaMP metabolomic pathway and metabolite database.
    """

    url = urls.urls['ramp']['url']
    c = curl.Curl(url, large = True, silent = False, compr = 'gz')

    return c._gzfile_mode_r


def ramp_raw(
        tables: list[str] = None,
        sqlite: bool = False,
        **kwargs
    ) -> dict[str, pd.DataFrame, sqlite3.Connection]:
    """
    Retrieve RaMP database contents from raw SQL dump.

    Args:
        tables:
            One or more tables to retrieve. If None, all tables are retrieved.
        sqlite:
            Return an SQLite database instead of a pandas DataFrame.
        kwargs:
            Options for the SQLite database: this way you can point to a new
            or existing database, while by default, an in-memory, temporary
            database is used.

    Returns:
        Either a dictionary with the table names as keys and  pandas dataframes
        as values, or an SQLite database connection.
    """

    fp = _ramp_sqldump()

    return sqldump.tables(
        sqldump = fp,
        tables = tables,
        return_df = True,
        return_sqlite = sqlite,
        con_param = kwargs,
        source_id = (fp.name, f'{os.path.getmtime(fp.name):.0f}'),
    )


def ramp_list_tables() -> dict[str, list[str]]:
    """
    List the tables of the RaMP database from SQL dump.
    """

    return sqldump.list_tables(_ramp_sqldump())


def ramp_show_tables() -> None:
    """
    Show the tables of the RaMP database from SQL dump.
    """

    pprint.pprint(ramp_list_tables())


def ramp_mapping(
        id_type_a: str,
        id_type_b: str,
        return_df: bool = False,
        curies: bool = False,
    ) -> dict[str, set[str]] | pd.DataFrame:
    """
    Retrieve the mapping between two identifiers.

    Args:
        id_type_a:
            The identifier type of the first identifier.
        id_type_b:
            The identifier type of the second identifier.
        return_df:
            Return a pandas DataFrame instead of a dictionary.
        curies:
            Do not remove CURIEs from the identifiers.

    Returns:
        A dictionary with the mapping between the two identifiers.
    """

    query = (
        'SELECT DISTINCT a.sourceId as id_type_a, b.sourceId as id_type_b '
        'FROM '
        '   (SELECT sourceId, rampId '
        '    FROM source '
        f'   WHERE geneOrCompound = "compound" AND IDtype = "{id_type_a}") a '
        'JOIN '
        '   (SELECT sourceId, rampId '
        '    FROM source '
        f'   WHERE geneOrCompound = "compound" AND IDtype = "{id_type_b}") b '
        'ON a.rampId = b.rampId;'
    )

    con = ramp_raw(tables = 'source', sqlite = True)
    df = pd.read_sql_query(query, con)

    if not curies:

        df[df.columns] = df[df.columns].apply(
            lambda y: [x.split(':', maxsplit = 1)[-1] for x in y],
        )

    return (
        df
            if return_df else
        df.groupby('id_type_a')['id_type_b'].apply(set).to_dict()
    )


def ramp_id_types(
        entity_type: Literal['gene', 'compound'] | None = None,
    ) -> set[str]:
    """
    List the identifier types of the RaMP database.
    """

    query = (
        'SELECT DISTINCT(s.IDtype) as id_type FROM source s' +
        (f' WHERE geneOrCompound = "{entity_type}";' if entity_type else ';')
    )
    con = ramp_raw(tables = 'source', sqlite = True)
    df = pd.read_sql_query(query, con)

    return set(df['id_type'])


def ramp_id_types_2(
        entity_type: Literal['gene', 'compound'] | None = None,
    ) -> set[str]:
    """
    List the identifier types of the RaMP database.

    Same output as `ramp_id_types`, but works by the API while the former
    extracts the data from the MySQL dump. The API means a fast, small
    download, while the SQL dump is huge and slow to process, but might
    be already available in the cache.
    """

    entity_types = {
        'compound': 'Metabolites',
        'gene': 'Genes/Proteins',
    }

    url = urls.urls['ramp']['api'] % 'id-types'
    c = curl.Curl(url, silent = True, large = False)

    return {
        id_type.strip()
        for i in json.loads(c.result)['data']
        if not entity_type or i['analyteType'] == entity_types[entity_type]
        for id_type in i['idTypes'].split(',')
    }
