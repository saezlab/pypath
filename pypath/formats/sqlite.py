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

"""
Utility functions for working with SQLite databases.
"""
from __future__ import annotations

import sqlite3
from pathlib import Path
from typing import Callable, TypeVar, BinaryIO

import pandas as pd
from pypath_common import _misc
from pypath.share import session, cache

__all__ = [
    'table_names',
    'list_columns',

    'sqlite_cache_path',
    'raw_tables',
    'iter',
]

_log = session.Logger(name = 'sqlite')._log


def table_names(con: sqlite3.Connection) -> list[str]:
    """
    Retrieves the names of all tables from a SQLite database.

    Args:
        con: An open SQLite database connection.

    Returns:
        A list of table names as strings.
    """
    cur = con.cursor()
    cur.execute('SELECT name FROM sqlite_master WHERE type = "table"')

    return [row[0] for row in cur.fetchall()]


def list_columns(con: sqlite3.Connection) -> dict[str, list[str]]:
    """
    Retrieves a list of tables and their columns from a SQLite database.

    Args:
        con: An open SQLite database connection.

    Returns:
        A dictionary mapping table names to a list of their column names.
    """
    return {
        table: [col[1] for col in con.execute(f'PRAGMA table_info({table})')]
        for table in table_names(con)
    }


def sqlite_cache_path(database: str, version: str) -> Path:
    """Constructs the standard local cache path for a SQLite database.

    The path is generated based on the pypath cache directory, the database
    name, and its version.

    Args:
        database: The name of the database (e.g., 'ChEMBL').
        version: The version string for the database.

    Returns:
        A `pathlib.Path` object for the cached SQLite file.
    """
    return Path(
        cache.get_cachedir()) / f'{database}_SQLite_{version}.sqlite'


DownloadResult = TypeVar("DownloadResult")


def download_sqlite(
        download_callback: Callable[[], DownloadResult],
        database: str,
        extractor: Callable[[DownloadResult], BinaryIO],
        version: str,
        connect: bool = True,
    ) -> sqlite3.Connection | Path:
    """
    Downloads and caches a database in SQLite format if not already present.

    This generic function checks for the existence of a specific version of a
    database in the local cache. If the file is not found, it uses the
    provided callback functions to download and extract it.

    Args:
        download_callback: A callable that performs the download and returns
            an object containing the downloaded data stream.
        database: The name of the database (e.g., 'RaMP').
        extractor: A callable that takes the result from `download_callback`
            and returns a readable file-like object for the SQLite database.
        version: The version of the database to download.
        connect: If True, returns an open `sqlite3.Connection` object.
            If False, returns the `pathlib.Path` to the database file.

    Returns:
        An open SQLite database connection or the path to the file.
    """
    sqlite_path = sqlite_cache_path(database, version)
    exists = Path.exists(sqlite_path)
    _log(f'{database} SQLite path: `{sqlite_path}`; exists: {exists}.')

    if not exists:

        # Resource specific download
        c = download_callback()

        # Resource specific sqlite file extractor
        file_stream = extractor(c)

        with open(sqlite_path, 'wb') as fp:

            # Reads 1 MB at a time
            while chunk := file_stream.read(1024 ** 2):

                fp.write(chunk)

    if connect:

        _log(f'Opening SQLite: `{sqlite_path}`')
        con = sqlite3.connect(sqlite_path)

        return con

    return sqlite_path


def _create_in_memory_subset_db(
        full_db_path: Path,
        tables: str |list[str],
        **kwargs,
    ) -> sqlite3.Connection:
    """
    Creates a new in-memory SQLite database and populates it with specified tables
    from a full database file.

    Args:
        full_db_path: Path to the full SQLite database file.
        tables: A list of table names to copy into the in-memory database.
        **kwargs: Additional keyword arguments to pass to `sqlite3.connect`
            for the in-memory database.

    Returns:
        An open `sqlite3.Connection` object to the new in-memory database.
    """
    # instructs sqlite3 to create a new in-memory database (RAM)
    _con_param = {'database': ':memory:'}

    _con_param.update(kwargs or {})

    # Creates the in-memory database and copies the requested tables
    con = sqlite3.connect(**_con_param)
    cur = con.cursor()
    cur.execute(f"ATTACH DATABASE '{full_db_path}' AS full")

    # Copy the requested tables
    for table in tables:
        cur.execute(f"CREATE TABLE {table} AS SELECT * FROM full.{table}")

    cur.execute('DETACH DATABASE full')
    con.commit()

    return con


def _extract_tables_to_data_structures(
        full_db_con: sqlite3.Connection,
        tables: str |list[str],
        return_df: bool,
    ) -> dict[str, list[tuple] | pd.DataFrame] | list[tuple] | pd.DataFrame:
    """
    Extracts data from specified tables in a SQLite database into DataFrames
    or lists of tuples.

    Args:
        full_db_con: An open SQLite database connection to the full database.
        tables: A list of table names to extract.
        return_df: If True, returns data as pandas DataFrames. Otherwise,
            returns lists of tuples.

    Returns:
        A dictionary mapping table names to DataFrames/lists of tuples, or a single
        DataFrame/list if only one table was requested.
    """

    # Default SQL query
    q = 'SELECT * FROM %s'

    if return_df:
        callback = lambda table_name: pd.read_sql_query(q % table_name, full_db_con)
    else:
        cur = full_db_con.cursor()
        callback = lambda table_name: cur.execute(q % table_name).fetchall()

    result = {t: callback(t) for t in tables}

    # If only one table was requested, return a single DataFrame/table
    if len(result) == 1:
        result = _misc.first(result.values())

    return result


def raw_tables(
        sqlite_callback: Callable[..., sqlite3.Connection | Path],
        tables: str | list[str] | None = None,
        sqlite: bool = False,
        return_df: bool = True,
        **kwargs,
    ) -> dict[str, list[tuple] | pd.DataFrame] | list[tuple] | pd.DataFrame | sqlite3.Connection:
    """
    Retrieves contents from a specified SQLite database.

    This function provides a flexible way to query one or more tables from a
    SQLite database, returning the data in various formats.

    Args:
        sqlite_callback: A callable that returns an open connection to the
            target SQLite database. Arguments are dependent on the specific
            database being accessed.
        tables: A string or list of strings specifying the table(s) to
            retrieve. If None, all tables are retrieved.
        sqlite: If True, returns a new in-memory SQLite database containing
            the requested tables.
        return_df: If True, returns data as pandas DataFrames. Otherwise,
            returns lists of tuples.
        **kwargs: Additional keyword arguments to pass to `sqlite3.connect`
            when creating a new in-memory database (if `sqlite=True`).

    Returns:
        A dictionary of table names to DataFrames/lists, or a connection.
    """
    full_db_path = sqlite_callback(connect=False)
    full_db_con = sqlite3.connect(full_db_path)

    if tables is None:
        tables = list(table_names(full_db_con))

    # Creates a SQLite database in memory with the requested tables
    if sqlite:
        return _create_in_memory_subset_db(full_db_path, tables, **kwargs)
    # Extracts the requested tables to data structures (DataFrames or lists)
    else:
        result = _extract_tables_to_data_structures(full_db_con, tables, return_df)
        full_db_con.close()
        return result


def iter(sqlite_callback: callable, table: str) -> Generator[tuple[Any]]:
    """
    Yields rows from a table in a SQLite database.

    This generator function is memory-efficient for iterating over large
    tables, as it retrieves one row at a time.

    Args:
        sqlite_callback: A callable that returns an open connection to the
            target SQLite database.
        table: The name of the table to iterate over.

    Yields:
        A tuple representing a single row from the table.
    """
    con = sqlite_callback()
    cur = con.cursor()

    cur.execute(f'SELECT * FROM {table}')

    for row in cur:

        yield row

    con.close()
