#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2023
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: Dénes Türei (turei.denes@gmail.com)
#           Nicolàs Palacio
#           Sebastian Lobentanzer
#           Erva Ulusoy
#           Olga Ivanova
#           Ahmet Rifaioglu
#           Sebastian Lobentanzer
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

from __future__ import annotations

"""
Access the RaMP metabolomic pathway and metabolite database.
"""

from typing import IO, Iterable
import os
import collections
import re

import sqlparse
# import pymysql
import pandas as pd

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.common as common
import pypath.share.session as session

_log = session.Logger(name = 'ramp_input')._log


def _ramp_sqldump() -> IO:
    """
    Download the RaMP metabolomic pathway and metabolite database.
    """

    url = urls.urls['ramp']['url']
    c = curl.Curl(url, large = True, silent = False, compr = 'gz')
    return c._gzfile_mode_r


def ramp_raw(tables: list[str] = None) -> dict[str, pd.DataFrame]:
    """
    Retrieve RaMP database contents from raw SQL dump.

    Args:
        tables:
            One or more tables to retrieve.

    Returns:
        A dictionary with the table names as keys and pandas dataframes as values.
    """

    return _sqldump_tables(_ramp_sqldump(), tables, return_df = True)


def ramp_list_tables() -> dict[str, list[str]]:
    """
    List the tables of the RaMP database from SQL dump.
    """

    return _sqldump_list_tables(_ramp_sqldump())


def _sqldump_tables(
        sqldump: str | IO,
        tables: str | Iterable[str] | None,
        return_df: bool = True,
        **kwargs
    ) -> tuple[dict, dict] | dict[str, pd.DataFrame]:
    """
    From a SQL dump, retrieve tables as data frames or as lists of tuples.

    Args:
        sqldump: SQL dump file. Path or file-like object.
        table: Name of the table.
        return_df: If True, return dict of pandas dataframes.
        **kwargs: Passed to `pypath.share.curl.FileOpener`.

    Returns:
        Contents of the table.
    """

    sqldump = _sqldump_open(sqldump, **kwargs)
    tables = common.to_set(tables) or _sqldump_list_tables(sqldump)

    contents = {}
    headers = {}

    for table in tables:

        _sqldump_seek(sqldump, table)
        the_table = _sqldump_one_table(sqldump, return_df)
        contents[table] = the_table if return_df else the_table[1]
        headers[table] = None if return_df else the_table[0]

        if not tables - set(contents.keys()):

            break

    return contents if return_df else headers, contents


def _sqldump_one_table(
        sqldump: IO,
        return_df: bool = True,
    ) -> tuple[list[str], list[tuple]] | pd.DataFrame:

    parser = sqlparse.parsestream(sqldump)
    header = []
    content = []
    inserts = 0

    for statement in parser:

        if statement.get_type() == 'INSERT':

            inserts += 1
            sublists = statement.get_sublists()
            table_info = next(sublists)
            table_name = table_info.get_name()

            if inserts == 1:
                
                the_table_name = table_name
                _log(f'Processing table: {table_name}')

            if table_name != the_table_name:

                break

            header = [
                col.get_name()
                for col in table_info.get_parameters()
            ]

            content.extend(
                tuple(
                    s.value.strip('"\'')
                    for s in next(rec.get_sublists()).get_identifiers()
                )
                for rec in next(sublists).get_sublists()
            )

    _log(
        f'Processed {len(content)} records in {inserts} INSERT statements '
        f'from table {table_name}.'
    )

    if return_df:

        df = pd.DataFrame.from_records(content, columns = header)

        _log(
            f'Data frame of table {table_name}: '
            f'{df.shape[0]} records, {common.df_memory_usage(df)}.'
        )

        return df

    else:

        return header, content


def _sqldump_table_2(
        sqldump: str | IO,
        tables: str | Iterable[str],
        return_df: bool = True,
        **kwargs
    ) -> list[tuple]:
    """
    From a SQL dump, retrieve a table as list of tuples.

    Args:
        sqldump: SQL dump file. Path or file-like object.
        table: Name of the table.
        return_df: If True, return dict of pandas dataframes.
        **kwargs: Passed to `pypath.share.curl.FileOpener`.

    Returns:
        Contents of the table.
    """

    tables = common.to_set(tables)
    sqldump = _sqldump_open(sqldump, **kwargs)

    con = pymysql.connect(
        
    )

    contents = collections.defaultdict(list)
    headers = {}





def _sqldump_list_tables(sqldump: str | IO) -> dict[str, list[str]]:
    """
    From a SQL dump, retrieve a list of table names.
    """
    
    recreate = re.compile(r'^CREATE TABLE `(.*)` \($')
    recol = re.compile(r'^\s*`(.*)`')
    tables = collections.defaultdict(list)
    current_table = None

    sqldump = _sqldump_open(sqldump)

    for line in sqldump:

        match = recreate.match(line)

        if match:

            current_table = match.group(1)

        match = recol.match(line)

        if match:

            tables[current_table].append(match.group(1))

        if line.strip() == ')':
            
            current_table = None

    return dict(tables)


def _sqldump_open(sqldump: str | IO, **kwargs) -> IO:
    """
    Open a SQL dump file.

    Args:
        sqldump: SQL dump file. Path or file-like object.

    Returns:
        The SQL dump file opened for reading, pointer at zero bytes.
    """

    if isinstance(sqldump, str) and os.path.exists(sqldump):

        fo = curl.FileOpener(sqldump, large = True, **kwargs)
        sqldump = fo.fileobj

    sqldump.seek(0)

    return sqldump


def _sqldump_seek(sqldump: IO, table: str, insert: bool = False) -> None:
    """
    Seek to a table in a SQL dump file.

    Args:
        sqldump:
            SQL dump file.
        table:
            Name of the table.
        insert:
            If True, seek to the beginning of the INSERT statement.
    """

    seek_for = f'INSERT INTO `{table}`' if insert else f'DROP TABLE IF EXISTS `{table}`'
    beginning_of_line = 0

    line = sqldump.readline()

    while line:

        if line.startswith(seek_for):

            sqldump.seek(beginning_of_line)
            return

        beginning_of_line = sqldump.tell()
        line = sqldump.readline()

    sqldump.seek(0)
    _sqldump_seek(sqldump, table, insert)
