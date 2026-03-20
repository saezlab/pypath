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

All public functions accept an ``sqlite`` parameter (default ``True``)
that selects between the SQLite build and the legacy SQL dump backend.
Functions that are only available via the SQLite backend do not have
this parameter.
"""

from typing import TYPE_CHECKING

if TYPE_CHECKING:

    import sqlite3

from collections.abc import Generator

import pandas as pd

from . import _sqldump
from . import _sqlite as _sqlite_mod
from ._sqlite import ramp_sqlite, ramp_omnipathmetabo
from ._rest import ramp_id_types as _rest_id_types
from ._mapping import ramp_mapping, ramp_synonym_mapping


def ramp_raw(
        tables: str | list[str] | None = None,
        return_sqlite: bool = False,
        return_df: bool = True,
        sqlite: bool = True,
        **kwargs,
    ) -> dict[str, list[tuple] | pd.DataFrame] | list[tuple] | pd.DataFrame | sqlite3.Connection:
    """
    Retrieve RaMP database contents.

    Args:
        tables:
            One or more tables to retrieve. If None, all tables are
            retrieved.
        return_sqlite:
            Return an SQLite database connection instead of data
            structures.
        return_df:
            Return pandas DataFrames. Ignored when ``return_sqlite``
            is True.
        sqlite:
            Use the SQLite build (default). If False, use the legacy
            SQL dump.
        kwargs:
            Extra options passed to the SQL dump backend's SQLite
            connection (ignored when ``sqlite=True``).

    Returns:
        Either a dictionary with table names as keys and pandas
        dataframes as values, or an SQLite database connection.
    """

    if sqlite:

        return _sqlite_mod.raw(
            tables = tables,
            return_sqlite = return_sqlite,
            return_df = return_df,
        )

    else:

        return _sqldump.ramp_raw(
            tables = tables,
            sqlite = return_sqlite,
            return_df = return_df,
            **kwargs,
        )


def ramp_list_tables(sqlite: bool = True) -> dict[str, list[str]]:
    """
    List the tables and their columns from the RaMP database.

    Args:
        sqlite:
            Use the SQLite build (default). If False, use the legacy
            SQL dump.

    Returns:
        A dictionary mapping table names to lists of column names.
    """

    if sqlite:

        return _sqlite_mod.list_tables()

    else:

        return _sqldump.ramp_list_tables()


def ramp_show_tables(sqlite: bool = True):
    """
    Pretty print the tables of the RaMP database.

    Args:
        sqlite:
            Use the SQLite build (default). If False, use the legacy
            SQL dump.
    """

    if sqlite:

        _sqlite_mod.show_tables()

    else:

        _sqldump.ramp_show_tables()


def ramp_table_fields(table: str) -> list[str]:
    """
    List the column names of a single table in the RaMP database.

    Args:
        table:
            Name of the table.

    Returns:
        A list of column names.
    """

    return _sqlite_mod.table_fields(table)


def ramp_id_types(
        entity_type: str | None = None,
    ) -> set[str]:
    """
    List the identifier types of the RaMP database.

    Uses the RaMP REST API to avoid downloading the full database.

    Args:
        entity_type:
            Filter by entity type, e.g. ``'gene'`` or ``'compound'``.

    Returns:
        A set of identifier type strings.
    """

    return _rest_id_types(entity_type = entity_type)


def ramp_table(table: str) -> Generator[tuple]:
    """
    Iterate over the rows of a table in the RaMP database.

    Args:
        table:
            Name of the table.

    Yields:
        Tuples representing individual rows.
    """

    return _sqlite_mod.table(table)


def ramp_table_df(table: str) -> pd.DataFrame:
    """
    Retrieve a table from the RaMP database as a data frame.

    Args:
        table:
            Name of the table.

    Returns:
        The full table as a pandas DataFrame.
    """

    return _sqlite_mod.table_df(table)
