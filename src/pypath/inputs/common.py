#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2020
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#                  Nicolàs Palacio
#                  Olga Ivanova
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

from past.builtins import xrange, range

import sys
import xlrd
from xlrd.biffh import XLRDError

import pypath.share.session as session_mod

_logger = session_mod.Logger(name = 'dataio')
_log = _logger._log
_console = _logger._console

if 'unicode' not in __builtins__: unicode = str


def read_xls(xls_file, sheet = '', csv_file = None, return_table = True):
    """
    Generic function to read MS Excel XLS file, and convert one sheet
    to CSV, or return as a list of lists
    """
    try:
        if hasattr(xls_file, 'read'):
            book = xlrd.open_workbook(
                file_contents = xls_file.read(),
                on_demand = True,
            )
        else:
            book = xlrd.open_workbook(xls_file, on_demand = True)
        try:
            sheet = book.sheet_by_name(sheet)
        except xlrd.biffh.XLRDError:
            sheet = book.sheet_by_index(0)
        table = [[unicode(c.value) for c in sheet.row(i)]
                 for i in xrange(sheet.nrows)]
        if csv_file:
            with open(csv_file, 'w') as csv:
                csv.write('\n'.join(['\t'.join(r) for r in table]))
        if not return_table:
            table = None
        book.release_resources()
        return table
    except IOError:
        sys.stdout.write('No such file: %s\n' % xls_file)
        sys.stdout.flush()


def csv_sep_change(csv, old, new):

    clean_csv = []
    bw_quotes = False

    for char in csv:
        if char == '\r':
            continue
        elif char == '"':
            bw_quotes = not bw_quotes
        elif char == '\n':
            if not bw_quotes:
                clean_csv.append(char)
            else:
                clean_csv.append(' ')
        elif char == old:
            if bw_quotes:
                clean_csv.append(char)
            else:
                clean_csv.append(new)
        else:
            clean_csv.append(char)

    return ''.join(clean_csv)


def _try_isoform(name):

    name = name.split('-')

    if len(name) > 1 and name[1].isdigit():

        isoform = int(name[1])
        main = name[0]

    else:

        main = '-'.join(name)
        isoform = None

    return main, isoform
