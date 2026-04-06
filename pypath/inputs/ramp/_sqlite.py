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
import re
import sqlite3

import pandas as pd

from pypath_common import _misc
import pypath.resources.urls as urls
from pypath.share.downloads import dm
import pypath.share.cache as cache
from pypath.formats import sqlite as _sqlite
from ._common import _log, _show_tables


_FALLBACK_VERSION = '2.5.4'

_INDEXES = (
    ('idx_source_rampId', 'source', 'rampId'),
    ('idx_source_IDtype_geneOrCompound', 'source', 'IDtype, geneOrCompound'),
    ('idx_analytesynonym_rampId', 'analytesynonym', 'rampId'),
    ('idx_metabolite_class_ramp_id', 'metabolite_class', 'ramp_id'),
    ('idx_chem_props_ramp_id', 'chem_props', 'ramp_id'),
)


def _ensure_indexes(con: sqlite3.Connection):
    """
    Create indexes on the RaMP database if they do not exist.
    """

    cur = con.cursor()

    for name, table, columns in _INDEXES:

        cur.execute(
            f'CREATE INDEX IF NOT EXISTS {name} ON {table} ({columns})'
        )

    con.commit()


def _latest_version() -> str:
    """
    Retrieve the latest RaMP database version from GitHub.

    Fetches the RaMP-DB GitHub repository db directory listing
    and extracts the highest version number from the SQLite file
    names.

    Returns:
        The latest version string, or the fallback version if
        the lookup fails.
    """

    try:

        import requests as _requests

        url = urls.urls['ramp']['github_db']
        resp = _requests.get(url, timeout=30)
        resp.raise_for_status()
        page_content = resp.text

        versions = sorted(
            re.findall(r'RaMP_SQLite_v([\d.]+)\.sqlite\.gz', page_content),
        )

        if versions:

            version = versions[-1]
            _log(f'Latest RaMP version from GitHub: {version}')

            return version

        _log('No RaMP versions found on GitHub, using fallback.')

    except Exception as e:

        _log(f'Failed to retrieve latest RaMP version: {e}')

    return _FALLBACK_VERSION


def ramp_sqlite(
        version: str | None = None,
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
            If None, the latest version is looked up from GitHub.

    Returns
        A SQLite database connection.
    """

    version = version or _latest_version()
    url = urls.urls['ramp']['github_sqlite'] % version

    import gzip
    import requests
    from pathlib import Path

    sqlite_path = _sqlite.sqlite_cache_path('RaMP', version)

    if not sqlite_path.exists():
        gz_path = Path(str(sqlite_path) + '.gz')

        if not gz_path.exists():
            _log(f'Downloading RaMP SQLite from {url}')
            gz_path.parent.mkdir(parents=True, exist_ok=True)
            resp = requests.get(url, stream=True, timeout=600)
            resp.raise_for_status()

            with open(gz_path, 'wb') as fp:
                for chunk in resp.iter_content(chunk_size=1024 ** 2):
                    fp.write(chunk)

        _log(f'Decompressing RaMP SQLite...')
        with gzip.open(gz_path, 'rb') as gz_in:
            with open(sqlite_path, 'wb') as out:
                while chunk := gz_in.read(1024 ** 2):
                    out.write(chunk)

        _log(f'Decompressed RaMP SQLite to {sqlite_path}')

    if connect:
        result = sqlite3.connect(sqlite_path)
    else:
        result = sqlite_path

    if connect:

        _ensure_indexes(result)

    return result


def raw(
        tables: str | list[str] | None = None,
        return_sqlite: bool = False,
        return_df: bool = True,
    ) -> dict[str, list[tuple] | pd.DataFrame] | list[tuple] | pd.DataFrame | sqlite3.Connection:
    """
    Retrieve RaMP database contents from its SQLite build.

    Args:
        tables:
            One or more tables to retrieve. If None, all tables are
            retrieved.
        return_sqlite:
            Return an SQLite database instead of a pandas DataFrame.
        return_df:
            Return a pandas data frame.

    Returns:
        Either a dictionary with the table names as keys and pandas
        dataframes as values, or an SQLite database connection.
    """

    return _sqlite.raw_tables(
        sqlite_callback = ramp_sqlite,
        tables = tables,
        sqlite = return_sqlite,
        return_df = return_df,
    )


def list_tables() -> dict[str, list[str]]:
    """
    List the tables and their columns from the RaMP SQLite database.

    Returns:
        A dictionary mapping table names to lists of column names.
    """

    con = ramp_sqlite()
    result = _sqlite.list_columns(con)
    con.close()

    return result


def show_tables():
    """
    Pretty print the tables of the RaMP SQLite database.
    """

    _show_tables(list_tables())


def table_fields(table: str) -> list[str]:
    """
    List the column names of a single table in the RaMP SQLite database.

    Args:
        table:
            Name of the table.

    Returns:
        A list of column names.
    """

    con = ramp_sqlite()
    columns = [
        col[1]
        for col in con.execute(f'PRAGMA table_info({table})')
    ]
    con.close()

    return columns


def id_types(
        entity_type: str | None = None,
    ) -> set[str]:
    """
    List the identifier types in the RaMP SQLite database.

    Args:
        entity_type:
            Filter by entity type, e.g. 'gene' or 'compound'.

    Returns:
        A set of identifier type strings.
    """

    query = (
        'SELECT DISTINCT(IDtype) AS id_type FROM source' +
        (f' WHERE geneOrCompound = "{entity_type}"' if entity_type else '')
    )
    con = ramp_sqlite()
    df = pd.read_sql_query(query, con)
    con.close()

    return set(df['id_type'])


def table(table: str) -> Generator[tuple]:
    """
    Iterate over the rows of a table in the RaMP SQLite database.

    Args:
        table:
            Name of the table.

    Yields:
        Tuples representing individual rows.
    """

    return _sqlite.iter(ramp_sqlite, table)


def table_df(table: str) -> pd.DataFrame:
    """
    Retrieve a table from the RaMP SQLite database as a data frame.

    Args:
        table:
            Name of the table.

    Returns:
        The full table as a pandas DataFrame.
    """

    con = ramp_sqlite()
    df = pd.read_sql_query(f'SELECT * FROM {table}', con)
    con.close()

    return df


def ramp_omnipathmetabo():
    """
    Retrieve RaMP database contents from its SQLite build for the
    omnipath.metabo database build. This function returns a namedtuple
    from a SQL query that joins the chemprops, metabolite_class,
    synonyms and source tables. These are the required field for the
    omnipath.metabo database build.

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
                    mw,
                    monoisotop_mass,
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
                        SELECT ramp_id,
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
                ORDER BY iso_smiles
                ''')

    col_names = [desc[0] for desc in cur.description]

    metabo_data = namedtuple('chem_props', col_names)

    for row in cur:

        yield metabo_data(*row)

    con.close()
