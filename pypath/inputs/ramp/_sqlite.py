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
from pypath.formats import sqlite as _sqlite
from ._common import _log, _show_tables


__all__ = [
    'ramp_sqlite',
    'ramp_sqlite_raw',
    'ramp_omnipathmetabo',
]

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

    url = urls.urls['ramp']['github_sqlite'] % version
    def _ramp_download() -> curl.Curl:
        """Callback to download the RaMP database."""
        return curl.Curl(urls,
                         large = True,
                         silent = False,
                         compr = 'gz'
                         )


    return _sqlite.download_sqlite(
        download_callback = _ramp_download,
        database = 'RaMP',
        version = version,
        connect = connect,
    )



def ramp_sqlite_raw(
        tables: str | list[str] | None = None,
        sqlite: bool = False,
        return_df: bool = True,
    ) -> dict[str, list[tuple] | pd.DataFrame] | list[tuple] | pd.DataFrame | sqlite3.Connection:
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

    return _sqlite.raw_tables(
        sqlite_callback = ramp_sqlite,
        tables = tables,
        sqlite = sqlite,
        return_df = return_df
    )


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