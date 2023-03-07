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

import sqlparse

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.common as common
import pypath.share.session as session

_log = session.Logger(name = 'ramp_input')._log


def ramp_raw(tables: list[str] = None) -> dict[str, list[tuple]]:
    """
    Retrieve raw SQL data.
    """

    url = urls.urls['ramp']['url']
    c = curl.Curl(url, large = True, silent = False, compr = 'gz')
    c.fileobj

    return c


def _sqldump_table(
        sqldump: str | IO,
        tables: str | Iterable[str],
        **kwargs
    ) -> list[tuple]:
    """
    From a SQL dump, retrieve a table as list of tuples.

    Args:
        sqldump: SQL dump file. Path or file-like object.
        table: Name of the table.
        **kwargs: Passed to `pypath.share.curl.FileOpener`.

    Returns:
        Contents of the table.
    """

    if isinstance(sqldump, str) and os.path.exists(sqldump):

        fo = curl.FileOpener(sqldump, large = True, **kwargs)
        sqldump = fo.fileobj

    tables = common.to_set(tables)
    sqldump.seek(0)

    contents = collections.defaultdict(list)
    headers = {}

    parser = sqlparse.parsestream(sqldump)

    for statement in parser:

        if not tables - set(contents.keys()):

            break

        if statement.get_type() == 'INSERT':

            sublists = statement.get_sublists()
            table_info = next(sublists)
            table_name = table_info.get_name()
            _log(f'Found table: {table_name}')

            if table_name not in tables:

                _log(f'Skipping table {table_name}.')
                continue
            
            _log(f'Processing table: {table_name}')

            headers[table_name] = [
                col.get_name()
                for col in table_info.get_parameters()
            ]

            contents[table_name].extend(
                tuple(
                    s.value.strip('"\'')
                    for s in next(rec.get_sublists()).get_identifiers()
                )
                for rec in next(sublists).get_sublists()
            )

    return headers, dict(contents)
