#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright 2014-2024
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
import pypath.share.common as common
import pypath.formats.sqldump as sqldump
from ._common import _log, _show_tables


def _ramp_sqldump() -> IO:
    """
    Download the RaMP metabolomic pathway and metabolite database.
    """

    url = urls.urls['ramp']['url']
    c = curl.Curl(url, large = True, silent = False, compr = 'gz')

    return c._gzfile_mode_r


def ramp_raw(
        tables: str | list[str] = None,
        sqlite: bool = False,
        return_df: bool = True,
        **kwargs
    ) -> dict[str, list[tuple] | pd.DataFrame | sqlite3.Connection]:
    """
    Retrieve RaMP database contents from raw SQL dump.

    Args:
        tables:
            One or more tables to retrieve. If None, all tables are retrieved.
        sqlite:
            Return an SQLite database instead of a pandas DataFrame.
        return_df:
            Return a pandas data frame.
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
        return_df = return_df,
        return_sqlite = sqlite,
        con_param = kwargs,
        source_id = (fp.name, f'{os.path.getmtime(fp.name):.0f}'),
    )


def ramp_list_tables() -> dict[str, list[str]]:
    """
    List the tables of the RaMP database from SQL dump.
    """

    return sqldump.list_tables(_ramp_sqldump())


def ramp_show_tables():
    """
    Show the tables of the RaMP database from SQL dump.
    """

    return _show_tables(ramp_list_tables())


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
