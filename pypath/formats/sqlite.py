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

import sqlite3

import pypath.share.session as session

_log = session.Logger(name = 'sqlite')._log


def list_tables(con: sqlite3.Connection) -> list[str]:
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
        for table in list_tables(con)
    }
