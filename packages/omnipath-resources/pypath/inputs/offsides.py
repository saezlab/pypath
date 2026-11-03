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

from typing import Generator

import collections

import pypath.share.curl as curl
import pypath.resources.urls as urls


def offsides_side_effects() -> Generator[tuple]:
    """
    Side effects from the OffSIDES (OFF label SIDE effectS) database.

    Retrieves individual drug side effect signals mined from the FDA's
    Adverse Event Reporting System.

    Yields:
        Tuples of drug side effect information.
    """

    return _sides_base(
        url_key = 'offsides',
        fields = (
            'drug_rxnorn',
            'drug',
            'condition_meddra',
            'condition',
            'prr',
            'prr_error',
            'mean_reporting_frequency',
        ),
        indices = (0, 1, 2, 3, 8, 9, 10),
        record_name = 'OffsideSideEffect',
    )


def _sides_base(
        url_key: str,
        fields: tuple[str],
        indices: tuple[int],
        record_name: str,
    ) -> Generator[tuple]:

    url = urls.urls[url_key]['url']
    c = curl.Curl(url, large = True, silent = False)

    result = set()
    record = collections.namedtuple(record_name, fields)

    _ = next(c.result)

    for line in c.result:

        line = line.strip().split(',')

        if not line:

            continue

        yield record(**{
            f: line[i].strip(' "')
            for f, i in zip(fields, indices)
        })
