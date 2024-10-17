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
Processing MySQL dump files.
"""

from typing import Generator, IO, Iterable
import os
import collections
import re
import sqlite3
import pickle

import sqlparse
import pandas as pd

import pypath.share.curl as curl
import pypath.share.common as common
import pypath.share.session as session
import pypath.share.cache as cache

_log = session.Logger(name = 'sqldump')._log


def tables(
        sqldump: str | IO,
        tables: str | Iterable[str] | None,
        return_df: bool = True,
        return_sqlite: bool = False,
        con_param: dict | None = None,
        source_id: str | tuple | None = None,
        **kwargs
    ) -> dict[str, list[tuple] | pd.DataFrame] | sqlite3.Connection:
    """
    From a SQL dump, retrieve tables as data frames or as lists of tuples.

    Args:
        sqldump:
            SQL dump file. Path or file-like object.
        table:
            Name of the table.
        return_df:
            If True, return dict of pandas dataframes.
        return_sqlite:
            If True, return dict of sqlite3.Connection objects.
            Has precedence over return_df.
        con_param:
            Connection parameters. By default we use an in-memory database,
            but you can also specify a path to a database file.
        source_id:
            Unique identifier for the source of the data. Necessary for
            caching, only if the path can not be determined from the
            file-like object.
        **kwargs: Passed to `pypath.share.curl.FileOpener`.

    Returns:
        Contents of the table.
    """

    common.log_memory_usage()
    source_id = get_source_id(sqldump, source_id)
    sqldump = open_sqldump(sqldump, **kwargs)
    tables = common.to_set(tables) or list_tables(sqldump)
    tables_done = set()

    if return_sqlite:

        _con_param = {'database': ':memory:'}
        _con_param.update(con_param or {})
        contents = sqlite3.connect(**_con_param)

    else:

        contents = {}

    for table in tables:

        content = table_from_cache(table, source_id)

        if not content:

            seek(sqldump, table)
            content = list(one_table(sqldump))
            table_to_cache(content, table, source_id)

            if not content:

                _log(f'Table `{table}` not found.')

        if return_df or return_sqlite:

            content = pd.DataFrame.from_records(
                content,
                columns = content[0]._fields,
            )

            _log(
                f'Data frame of table `{table}`: '
                f'{content.shape[0]} records, {common.df_memory_usage(content)}.'
            )
            common.log_memory_usage()

        if return_sqlite:

            _log(f'Loading table `{table}` into SQLite database.')
            content.to_sql(
                table,
                contents,
                if_exists = 'replace',
                index = False,
            )
            common.log_memory_usage()

        else:

            contents[table] = content

        tables_done.add(table)

        if not tables - tables_done:

            break

    if not return_sqlite and len(contents) == 1:

        contents = common.first(contents.values())

    return contents


def table_from_cache(
        table: str,
        source_id: tuple,
    ) -> list[tuple]:
    """
    Retrieve table contents from cache.
    """

    key = source_id + (table,)
    cache_path = cache.cache_path(key)

    if os.path.exists(cache_path):

        with open(cache_path, 'rb') as fp:

            _log(
                f'Loading table `{table}` from cache. '
                f'Cache path: `{cache_path}`.'
            )

            name, fields, content = pickle.load(fp)

            record = collections.namedtuple(name, fields)
            content = [record(*r) for r in content]

            return content

    return []


def table_to_cache(
        content: list[tuple],
        table: str,
        source_id: tuple,
    ) -> None:
    """
    Save table contents to cache.
    """

    key = source_id + (table,)
    cache_path = cache.cache_path(key)

    with open(cache_path, 'wb') as fp:

        _log(f'Saving table `{table}` to cache. Cache path: `{cache_path}`.')
        pickle.dump(
            (
                content[0].__class__.__name__,
                content[0]._fields,
                [tuple(r) for r in content]
            ),
            fp,
        )


def one_table(sqldump: IO) -> Generator[tuple]:
    """
    Retrieve one table from a MySQL dump.
    """

    parser = sqlparse.parsestream(sqldump)
    n_inserts = 0
    n_records = 0
    record = None
    table_name = the_table_name = None

    for statement in parser:

        if statement.get_type() == 'INSERT':

            n_inserts += 1
            sublists = statement.get_sublists()
            table_info = next(sublists)
            table_name = table_info.get_name()

            if n_inserts == 1:

                the_table_name = table_name
                _log(f'Processing table: `{table_name}`')

            if table_name != the_table_name:

                break

            if record is None:

                record = collections.namedtuple(
                    f'Ramp{table_name.capitalize()}',
                    [col.get_name() for col in table_info.get_parameters()],
                )
                globals()[record.__class__.__name__] = record

            for rec in next(sublists).get_sublists():

                n_records += 1

                yield record(*(
                    s.value.strip('"\'')
                    for s in next(rec.get_sublists()).get_identifiers()
                ))

    if table_name:

        _log(
            f'Processed {n_records} records in {n_inserts} INSERT statements '
            f'from table `{the_table_name}`.'
        )


def get_source_id(
        sqldump: str | IO,
        source_id: str | tuple | None = None,
    ) -> tuple:
    """
    Attempts to use the file name as source id.
    """

    source_id = (
        source_id or
        (
            os.path.basename(sqldump)
                if isinstance(sqldump, str) else
            getattr(sqldump, 'name', None)
        )
    )

    return common.to_tuple(source_id)


def endof_statement(line: str) -> bool:
    """
    Check if the current line is the end of a SQL statement.
    """

    return line.strip().endswith(';')


def list_tables(sqldump: str | IO) -> dict[str, list[str]]:
    """
    From a SQL dump, retrieve a list of table names.
    """

    recreate = re.compile(r'^CREATE TABLE `(.*)` \($')
    recol = re.compile(r'^\s*`(.*)`')
    tables = collections.defaultdict(list)
    current_table = None

    sqldump = open_sqldump(sqldump)

    for line in sqldump:

        line = line.strip()
        match = recreate.match(line)

        if match:

            current_table = match.group(1)

        match = recol.match(line)

        if match:

            tables[current_table].append(match.group(1))

        if endof_statement(line):

            current_table = None

    return dict(tables)


def open_sqldump(sqldump: str | IO, **kwargs) -> IO:
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


def seek(sqldump: IO, table: str, insert: bool = False) -> None:
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

    start_from = sqldump.tell()
    seek_for = (
        f'INSERT INTO `{table}`'
            if insert else
        f'DROP TABLE IF EXISTS `{table}`'
    )
    beginning_of_line = 0

    line = sqldump.readline()

    while line:

        if line.startswith(seek_for):

            sqldump.seek(beginning_of_line)
            return

        beginning_of_line = sqldump.tell()
        line = sqldump.readline()

    if start_from != 0:

        sqldump.seek(0)
        seek(sqldump, table, insert)

    else:

        _log(f'Table `{table}` not found in SQL dump.')
