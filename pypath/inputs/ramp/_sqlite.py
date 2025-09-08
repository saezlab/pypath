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

from collections.abc import Generator
from collections import namedtuple
from typing import Any
import os
import sqlite3

import pandas as pd

from pypath_common import _misc
import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.cache as cache
import pypath.formats.sqlite as _sqlite
from ._common import _log, _show_tables


__all__ = [
    'ramp_iter',
    'ramp_sqlite',
    'ramp_show_tables',
    'ramp_list_tables',
    'ramp_raw',
    'ramp_omnipathmetabo',
]


def _ramp_sqlite_path(version: str = '2.5.4') -> str:

    return os.path.join(
        cache.get_cachedir(),
        f'RaMP_SQLite_{version}.sqlite',
    )


def ramp_sqlite(
        version: str = '2.5.4',
        connect: bool = True,
    ) -> sqlite3.Connection | str:
    """
    Download the RaMP database in SQLite format.

    SQLite format builds are avialable in RaMP's GitHub repo:
        https://github.com/ncats/RaMP-DB/tree/main/db
    This function downloads and opens one of them.

    Args
        version:
            The version of the RaMP database to download.

    Returns
        A SQLite database connection.
    """

    sqlite_path = _ramp_sqlite_path(version)
    exists = os.path.exists(sqlite_path)
    _log(f'RaMP SQLite path: `{sqlite_path}`; exists: {exists}.')

    if not exists:

        url = urls.urls['ramp']['github_sqlite'] % version
        c = curl.Curl(url, large = True, silent = False, compr = 'gz')

        with open(sqlite_path, 'wb') as fp:

            while chunk := c.gzfile.read(1024 ** 2):

                fp.write(chunk)

    if connect:

        _log(f'Opening SQLite: `{sqlite_path}`')
        con = sqlite3.connect(sqlite_path)

        return con

    return sqlite_path


def ramp_list_tables() -> dict[str, list[str]]:
    """
    List the tables of the RaMP SQLite database.
    """

    con = ramp_sqlite()
    tables = _sqlite.list_columns(con)
    con.close()

    return tables


def ramp_show_tables():
    """
    Show the tables of the RaMP database from SQL dump.
    """

    return _show_tables(ramp_list_tables())


def ramp_raw(
        tables: str | list[str] = None,
        sqlite: bool = False,
        return_df: bool = True,
        prefixes: bool = True,
        **kwargs,
    ) -> dict[str, list[tuple] | pd.DataFrame | sqlite3.Connection]:
    """
    Retrieve RaMP database contents from its SQLite build.

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
        Either a dictionary with the table names as keys and pandas dataframes
        as values, or an SQLite database connection.
    """

    ramp_sqlite_path = ramp_sqlite(connect = False)
    tables = _misc.to_list(tables) or list(ramp_list_tables().keys())

    if sqlite:

        _con_param = {'database': ':memory:'}
        _con_param.update(kwargs or {})
        con = sqlite3.connect(**_con_param)
        cur = con.cursor()
        cur.execute(f"ATTACH DATABASE '{ramp_sqlite_path}' AS ramp")

        for table in tables:

            cur.execute(f"CREATE TABLE {table} AS SELECT * FROM ramp.{table}")

        cur.execute('DETACH DATABASE ramp')
        con.commit()

        return con

    else:

        con = sqlite3.connect(ramp_sqlite_path)
        q = f'SELECT * FROM %s'

        if return_df:

            callback = lambda x: pd.read_sql_query(q % x, con)

        else:

            cur = con.cursor()
            callback = lambda x: cur.execute(q % x).fetchall()

        result = {t: callback(t) for t in tables}

        if not prefixes and 'source' in result:
            result['source'].sourceId = result['source'].sourceId.apply(lambda x: x.split(':', maxsplit = 1)[1])

        if len(result) == 1:

            result = _misc.first(result.values())

        return result


def ramp_iter(table: str) -> Generator[tuple[Any]]:
    """
    Retrieve RaMP database contents from its SQLite build.

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
        Either a dictionary with the table names as keys and pandas dataframes
        as values, or an SQLite database connection.
    """

    con = ramp_sqlite()
    cur = con.cursor()

    cur.execute(f'SELECT * FROM {table}')

    for row in cur:

        yield row

    con.close()


def ramp_omnipathmetabo():
    """
    Retrieve RaMP database contents from its SQLite build for the omnipath.metabo database build.
    This function returns a namedtuple from a SQL query that joins the chemprops, metabolite_class,
    synonyms and source tables. These are the required field for the omnipath.metabo database build.

    Returns:
        A namedtuple out put of the SQL query.
    """

    con = ramp_sqlite()
    cur = con.cursor()

    cur.execute('''
                SELECT
                    chem_props.ramp_id AS ramp_id, 
                    chem_data_source,
                    chem_source_id,
                    iso_smiles,
                    inchi_key,
                    inchi,
                    CAST(mw AS DOUBLE) AS mw,
                    CAST(monoisotop_mass AS DOUBLE) AS monoisotop_mass,
                    common_name,
                    mol_formula,
                    synonyms,
                    classes,
                    sources
                FROM chem_props

                -- sub query to group the synonyms from the analytesynonym table
                LEFT JOIN (
                        SELECT rampId,
                               GROUP_CONCAT(synonym, ', ') AS synonyms
                        FROM analytesynonym
                        GROUP BY rampId
                        ) syn
                    ON chem_props.ramp_id = syn.rampId
                -- sub query to group the classes from the metabolite_class table
                LEFT JOIN (
                        SELECT DISTINCT ramp_id,
                               GROUP_CONCAT(class_source_id || '|' || class_level_name || '|' || class_name, ', ') AS classes
                        FROM metabolite_class
                        GROUP BY ramp_id
                        ) mc
                    ON chem_props.ramp_id = mc.ramp_id
                -- sub query to group the sources from the source table
                LEFT JOIN (
                        SELECT rampId,
                               GROUP_CONCAT(sourceID, ', ') AS sources
                        FROM source
                        GROUP BY rampId
                        ) src
                    ON chem_props.ramp_id = src.rampId
                
                GROUP BY 
                    chem_props.ramp_id,
                    chem_source_id
                ORDER BY iso_smiles
                ''')
    
    col_names = [desc[0] for desc in cur.description]

    metabo_data = namedtuple('chem_props', col_names)

    for row in cur:

        yield metabo_data(*row)

    con.close()

    