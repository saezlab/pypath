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

import re

import pandas as pd

import pypath.share.curl as curl
import pypath.share.common as common
import pypath.resources.urls as urls


def expasy_enzyme_classes(
        return_df: bool = False,
    ) -> list[tuple] | pd.DataFrame:
    """
    Enzyme classification from the ExPASy database.

    Args:
        return_df:
            Return a pandas data frame.

    Returns:
        Tuples with the first 3 digits and the names of EC classes.
    """

    url = urls.urls['expasy']['enzclass']

    c = curl.Curl(
        url,
        silent=False,
        large=True,
        encoding='utf-8',
        default_mode='r',
    )
    reec = re.compile(
        r'^(\d+)\.\s?'
        r'(?:(\d+))?-?\.\s*'
        r'(?:(\d+))?-?\.\s*'
        r'-?\s+'
        r'([-\w\(\)\s\+]+)\.$'
    )

    result = []

    for line in c.result:


        m = reec.match(line)

        if m:

            result.append(m.groups())

    return pd.DataFrame(result) if return_df else result


def expasy_enzymes(return_df: bool = False) -> dict[key, tuple] | pd.DataFrame:
    """
    Enzyme data from the ExPASy database.

    Args:
        return_df:
            Return a data frame.
    """

    url = urls.urls['expasy']['enzymes']

    c = curl.Curl(
        url,
        silent = False,
        large = True,
        encoding = 'utf-8',
        default_mode = 'r',
        slow = True,
    )

    reid = re.compile(r'(\d+\.\d+\.\d+\.[\dn]+)')
    reup = re.compile(r'([\w]+), ([\w]+_\w+)')
    reeq = re.compile(r'\(\d\)')
    recc = re.compile(r'-!-')

    result = []
    new_enzyme = lambda: {
        k: []
        for k in ('uniprots', 'entries', 'ca', 'cc', 'an')
    }
    enzyme = new_enzyme()

    for line in c.result:

        prefix, *line = line.split(maxsplit = 1)

        if prefix == 'ID':

            enzyme['ec'] = reid.match(line[0]).group(0)

        elif prefix == 'DE':

            enzyme['de'] = line[0].strip('.\n ')

        elif prefix in {'CA', 'CC', 'AN'} and line:

            enzyme[prefix.lower()].append(line[0].strip('. \n'))

        elif prefix == 'DR':

            uniprots, entries = list(zip(*reup.findall(line[0])))
            enzyme['uniprots'].extend(uniprots)
            enzyme['entries'].extend(entries)

        elif prefix == '//' and enzyme.get('ec', None):

            for key, rexp in zip(('ca', 'cc'), (reeq, recc)):

                enzyme[key] = common.del_empty(
                    x.strip()
                    for x in rexp.split(' '.join(enzyme[key]))
                )

            result.append(enzyme)
            enzyme = new_enzyme()

    return pd.DataFrame(result) if return_df else {e['ec']: e for e in result}
