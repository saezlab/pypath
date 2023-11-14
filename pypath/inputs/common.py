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

from future.utils import iteritems
from past.builtins import xrange, range

import os
import sys
import warnings
import json

from typing import Any, Callable, Dict, IO, List, Optional, Union

import xlrd
import openpyxl
import glom

import pypath.share.session as session_mod
import pypath.share.common as common
import pypath_common._constants as _const

_logger = session_mod.Logger(name = 'inputs_common')
_log = _logger._log
_console = _logger._console

if 'unicode' not in __builtins__: unicode = str


def read_xls(
        xls_file,
        sheet = 0,
        use_openpyxl = False,
        cell_range = None,
    ):
    """
    Generic function to read MS Excel XLS file, and convert one sheet
    to CSV, or return as a list of lists
    """

    table = []
    opened_here = False

    if isinstance(xls_file, str):

        if os.path.exists(xls_file):

            xls_file = open(xls_file, 'rb')
            opened_here = True

        else:

            raise FileNotFoundError(xls_file)

    if not use_openpyxl:

        try:

            _log('Reading XLS(X) by xlrd.')

            if hasattr(xls_file, 'read'):

                book = xlrd.open_workbook(
                    file_contents = xls_file.read(),
                    on_demand = True,
                )

            try:
                if isinstance(sheet, int):
                    sheet = book.sheet_by_index(sheet)
                else:
                    sheet = book.sheet_by_name(sheet)
            except xlrd.biffh.XLRDError:
                sheet = book.sheet_by_index(0)

            table = [
                [str(c.value) for c in sheet.row(i)]
                for i in xrange(sheet.nrows)
            ]

            use_openpyxl = False

        except IOError:

            raise FileNotFoundError(xls_file)

        except Exception as e:

            _log('Failed to read by xlrd, falling back to openpyxl.')
            _logger._log_traceback()
            use_openpyxl = True

    if use_openpyxl:

        try:

            _log('Reading XLS(X) by openpyxl.')

            book = openpyxl.load_workbook(
                filename = xls_file,
                read_only = True,
                data_only = True,
            )

        except Exception as e:

            _log(f'Failed to read `{xls_file}` by openpyxl.')
            _logger._log_traceback()
            raise ValueError('Could not open xls: %s' % xls_file)

        try:

            if type(sheet) is int:
                sheet = book.worksheets[sheet]
            else:
                sheet = book[sheet]

        except:

            sheet = book.worksheets[0]

        # this is to suppress the openpyxl unknown extension warnings
        # which we can not avoid as the xlsx files were produced not by us
        with warnings.catch_warnings():

            warnings.simplefilter('ignore')

            table = [
                [
                    (
                        cell
                            if isinstance(cell, str) else
                        cell.value
                            if cell is not None else
                        ''
                    )
                    for cell in row
                ]
                for row in (sheet[cell_range] if cell_range else sheet.values)
            ]

    if 'book' in locals() and hasattr(book, 'release_resources'):

        book.release_resources()

    if opened_here:

        xls_file.close()

    return table


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


def read_table(
        cols,
        fileObject = None,
        data = None,
        sep = '\t',
        sep2 = None,
        rem = None,
        hdr = None,
        encoding = 'ascii',
    ):
    """
    Generic function to read data tables.

    fileObject : file-like
        Any file like object: file opened for read, or StringIO buffer
    cols : dict
        Dictionary of columns to read. Keys identifying fields are returned
        in the result. Values are column numbers.
    sep : str
        Field separator of the file.
    sep2 : dict
        Subfield separators and prefixes.
        E.g. {2: ',', 3: '|'}
    hdr : int
        Number of header lines. If None, no headers assumed.
    rem : list
        Strings to remove. For each line these elements will be replaced with ''.
    """

    rem = rem or []

    if data is None:

        if hasattr(fileObject, 'readline'):

            fileObject.seek(0)

        if hdr:

            for h in xrange(0, hdr):

                _ = next(fileObject)

        data = fileObject

    else:

        data = [l.strip() for l in data.split('\n') if len(l) > 0][hdr:]

    res = []

    for l in data:

        if type(l) is bytes:

            l = l.decode(encoding)

        for r in rem:

            l = l.replace(r, '')

        l = [f.strip() for f in l.split(sep)]

        if len(l) > max(cols.values()):

            dic = {}

            for name, col in iteritems(cols):

                field = l[col].strip()

                _sep2 = (
                    sep2[col]
                        if isinstance(sep2, dict) and col in sep2 else
                    sep2
                        if isinstance(sep2, str) else
                    None
                )

                if _sep2:

                    field = tuple(
                        sf.strip()
                        for sf in field.split(_sep2)
                        if sf
                    )

                dic[name] = field

            res.append(dic)

    if fileObject is not None:

        fileObject.close()

    return res


def json_extract(
        data: Union[dict, list, str, IO],
        spec: dict,
    ) -> List[dict]:
    """
    Extracts fields of arbitrary depth from JSON data into a list of dicts.

    Args
        data: JSON as a string or a file-like object.
        spec: Dict of glom field specifications.
    """

    data = json_read(data)

    if isinstance(data, dict):

        data = [data]

    if not isinstance(data, list):

        msg = 'Don\'t know how to process data of type `%s`.' % type(data)
        raise TypeError(msg)


    return [
        glom.glom(rec, spec, default = _cons.GLOM_ERROR)
        for rec in data
    ]


def json_read(data: Union[str, IO, Any]) -> Union[list, dict, Any]:
    """
    Reads JSON from file or string, pass through for any other value.
    """

    if isinstance(data, IO):

        data = json.load(data)

    elif isinstance(data, str):

        data = json.loads(data)

    return data


GlomSpec = Union[str, tuple, dict, Callable]

GlomFields = Union[
    List[str],
    Dict[str, GlomSpec]
]

def glom_fields(fields: Optional[GlomFields] = None) -> Dict[str, GlomSpec]:
    """
    Generates a glom spec dict from a list or dict, protecting each field
    by glom.Coalesce.
    """

    fields = fields or {}

    fields = fields if isinstance(fields, dict) else dict(zip(fields, fields))

    fields = dict(
        (
            k,
            glom.Coalesce(v, default = None)
        )
        for k, v in fields.items()
    )

    return fields
