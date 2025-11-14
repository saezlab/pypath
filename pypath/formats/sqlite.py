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
Utility functions for working with SQLite databases.
"""

import os
import sqlite3

import pypath.share.session as session

__all__ = [
    'table_names',
    'list_columns',
    'list_tables',
    'sqlite_cache_path',
    'raw_tables',
    'iter',
]

_log = session.Logger(name = 'sqlite')._log


def table_names(con: sqlite3.Connection) -> list[str]:
    """
    From a SQLite database, retrieve a list of table names.
    """

    cur = con.cursor()
    cur.execute('SELECT name FROM sqlite_master WHERE type = "table"')

    return [row[0] for row in cur.fetchall()]


def list_columns(con: sqlite3.Connection) -> dict[str, list[str]]:
    """
    From a SQLite database, retrieve a list of column names.
    """

    return {
        table: [col[1] for col in con.execute(f'PRAGMA table_info({table})')]
        for table in table_names(con)
    }


def list_tables(con: sqlite3.Connection) -> dict[str, list[str]]:
    """
    List the tables of an SQLite database.
    """

    tables = list_columns(con)
    con.close()

    return tables


def sqlite_cache_path(database: str, version: str) -> str:

    return os.path.join(
        cache.get_cachedir(),
        f'{database}_SQLite_{version}.sqlite',
    )


def download_sqlite(
        download_callback: callable,
        database: str,
        version: str,
        connect: bool = True,
    ) -> sqlite3.Connection | str:
    """
    Download a database in SQLite format.

    Args
        version:
            The version of the database to download.

    Returns
        A SQLite database connection.
    """

    sqlite_path = sqlite_path(database, version)
    exists = os.path.exists(sqlite_path)
    _log(f'{database} SQLite path: `{sqlite_path}`; exists: {exists}.')

    if not exists:

        c = download_callback()

        with open(sqlite_path, 'wb') as fp:

            while chunk := c.gzfile.read(1024 ** 2):

                fp.write(chunk)

    if connect:

        _log(f'Opening SQLite: `{sqlite_path}`')
        con = sqlite3.connect(sqlite_path)

        return con

    return sqlite_path


def show_tables(con: sqlite3.Connection):
    """
    Show the tables of the RaMP database from SQL dump.
    """

    return pprint.pprint(list_tables(con))


def raw_tables(
        sqlite_callback: callable,
        tables: str | list[str] = None,
        sqlite: bool = False,
        return_df: bool = True,
        prefixes: bool = True,
        **kwargs,
    ) -> dict[str, list[tuple] | pd.DataFrame | sqlite3.Connection]:
    """
    Retrieve contents from an SQLite database.

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

    full_db_path = sqlite_callback(connect = False)
    full_db_con = sqlite_callback()
    tables = _misc.to_list(tables) or list(list_names(full_db_con))

    if sqlite:

        _con_param = {'database': ':memory:'}
        _con_param.update(kwargs or {})
        con = sqlite3.connect(**_con_param)
        cur = con.cursor()
        cur.execute(f"ATTACH DATABASE '{full_db_path}' AS full")

        for table in tables:

            cur.execute(f"CREATE TABLE {table} AS SELECT * FROM full.{table}")

        cur.execute('DETACH DATABASE full')
        con.commit()

        return con

    else:

        q = f'SELECT * FROM %s'

        if return_df:

            callback = lambda x: pd.read_sql_query(q % x, full_db_con)

        else:

            cur = full_db_con.cursor()
            callback = lambda x: cur.execute(q % x).fetchall()

        result = {t: callback(t) for t in tables}

        if not prefixes and 'source' in result:

            result['source'].sourceId = (
                result['source'].sourceId.apply(
                    lambda x: x.split(':', maxsplit = 1)[1]
                )
            )

        if len(result) == 1:

            result = _misc.first(result.values())

        return result


def iter(sqlite_callback: callable, table: str) -> Generator[tuple[Any]]:
    """
    Retrieve contents from an SQLite database.

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

    con = sqlite_callback()
    cur = con.cursor()

    cur.execute(f'SELECT * FROM {table}')

    for row in cur:

        yield row

    con.close()
